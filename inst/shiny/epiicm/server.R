##
## Server File for epiicm Shiny Application
##
## Run local: epiweb(class = "icm")
##

library(shiny)
library(EpiModel)

shinyServer(function(input, output, session) {

  # =========================================================================
  # Presets
  # =========================================================================
  # Track which modtype a preset expects, so we can distinguish
  # preset-driven modtype changes from manual ones
  preset_expected_modtype <- reactiveVal(NULL)

  observeEvent(input$preset, {
    if (input$preset == "Flu-like (SIR)") {
      preset_expected_modtype("SIR")
      updateSelectInput(session, "modtype", selected = "SIR")
      updateSliderInput(session, "inf.prob", value = 0.03)
      updateSliderInput(session, "act.rate", value = 10)
      updateSliderInput(session, "rec.rate", value = round(1 / 7, 3))
      updateNumericInput(session, "s.num", value = 500)
      updateNumericInput(session, "i.num", value = 1)
      updateNumericInput(session, "r.num", value = 0)
      updateNumericInput(session, "nsteps", value = 200)
      updateNumericInput(session, "nsims", value = 10)
      updateCheckboxInput(session, "enable_vital", value = FALSE)
      updateCheckboxInput(session, "enable_intervention", value = FALSE)
    } else if (input$preset == "STI-like (SIS)") {
      preset_expected_modtype("SIS")
      updateSelectInput(session, "modtype", selected = "SIS")
      updateSliderInput(session, "inf.prob", value = 0.2)
      updateSliderInput(session, "act.rate", value = 0.5)
      updateSliderInput(session, "rec.rate", value = 0.01)
      updateNumericInput(session, "s.num", value = 500)
      updateNumericInput(session, "i.num", value = 10)
      updateNumericInput(session, "nsteps", value = 500)
      updateNumericInput(session, "nsims", value = 10)
      updateCheckboxInput(session, "enable_vital", value = FALSE)
      updateCheckboxInput(session, "enable_intervention", value = FALSE)
    } else if (input$preset == "Measles-like (SIR)") {
      preset_expected_modtype("SIR")
      updateSelectInput(session, "modtype", selected = "SIR")
      updateSliderInput(session, "inf.prob", value = 0.5)
      updateSliderInput(session, "act.rate", value = 3)
      updateSliderInput(session, "rec.rate", value = 0.1)
      updateNumericInput(session, "s.num", value = 500)
      updateNumericInput(session, "i.num", value = 1)
      updateNumericInput(session, "r.num", value = 0)
      updateNumericInput(session, "nsteps", value = 100)
      updateNumericInput(session, "nsims", value = 10)
      updateCheckboxInput(session, "enable_vital", value = FALSE)
      updateCheckboxInput(session, "enable_intervention", value = FALSE)
    }
  })

  # Reset preset to "Custom" when user manually changes disease type.
  # If the modtype change matches what a preset just requested, consume
  # the expectation and leave the preset alone.
  observeEvent(input$modtype, {
    expected <- preset_expected_modtype()
    if (!is.null(expected) && input$modtype == expected) {
      preset_expected_modtype(NULL)
    } else if (input$preset != "Custom") {
      preset_expected_modtype(NULL)
      updateSelectInput(session, "preset", selected = "Custom")
    }
  })

  # =========================================================================
  # Dynamic UI: update plot type choices based on model type and vital
  # =========================================================================
  observe({
    type <- input$modtype
    vital <- isTRUE(input$enable_vital)
    choices <- c("Compartment Prevalence", "Compartment Size",
                 "Disease Incidence")
    if (type == "SIR") {
      choices <- c(choices, "Recovery Flow")
    }
    if (type == "SIS") {
      choices <- c(choices, "Re-susceptibility Flow")
    }
    if (vital) {
      choices <- c(choices, "Arrival Flow", "Departure Flows")
    }
    updateSelectInput(session, "compsel", choices = choices)
  })

  # Update death rate label based on differential checkbox
  observeEvent(input$diff_death_rates, {
    lbl <- if (isTRUE(input$diff_death_rates)) {
      "Death Rate, Susceptible"
    } else {
      "Death Rate"
    }
    updateNumericInput(session, "ds.rate", label = lbl)
  })


  # =========================================================================
  # Model state: reactiveValues to support both auto-run and button-run
  # =========================================================================
  rv <- reactiveValues(mod = NULL)

  # Shared function to build and run the model from current inputs
  build_and_run <- function() {
    s_num     <- input$s.num
    i_num     <- input$i.num
    r_num     <- input$r.num
    nsteps    <- input$nsteps
    nsims     <- input$nsims
    inf_prob  <- input$inf.prob
    act_rate  <- input$act.rate
    rec_rate  <- input$rec.rate
    modtype   <- input$modtype

    # Build param
    p_args <- list(inf.prob = inf_prob, act.rate = act_rate)
    if (modtype != "SI") {
      p_args$rec.rate <- rec_rate
    }

    # Intervention
    if (isTRUE(input$enable_intervention)) {
      p_args$inter.eff <- input$inter.eff
      p_args$inter.start <- as.numeric(input$inter.start)
    }

    # Vital dynamics
    if (isTRUE(input$enable_vital)) {
      p_args$a.rate <- input$a.rate
      p_args$ds.rate <- input$ds.rate
      if (isTRUE(input$diff_death_rates)) {
        p_args$di.rate <- input$di.rate
        if (modtype == "SIR") {
          p_args$dr.rate <- input$dr.rate
        }
      } else {
        p_args$di.rate <- input$ds.rate
        if (modtype == "SIR") {
          p_args$dr.rate <- input$ds.rate
        }
      }
    }

    param <- do.call(param.icm, p_args)

    # Build init — as.numeric() is critical: Shiny numericInput returns
    # integer class, but init.icm internally filters with
    # sapply(init, class) == "numeric" which excludes integers
    i_args <- list(s.num = as.numeric(s_num), i.num = as.numeric(i_num))
    if (modtype == "SIR") {
      i_args$r.num <- as.numeric(r_num)
    }
    init <- do.call(init.icm, i_args)

    # Build control
    control <- control.icm(type = modtype, nsteps = nsteps,
                            nsims = nsims, verbose = FALSE)

    # Run model
    tryCatch({
      icm(param, init, control)
    }, error = function(e) {
      showNotification(paste("Model error:", e$message),
                       type = "error", duration = 10)
      NULL
    })
  }

  # Run on button click
  observeEvent(input$runMod, {
    withProgress(message = "Running simulations...", {
      rv$mod <- build_and_run()
    })
  })

  # Auto-run once on startup after inputs are initialized
  observe({
    req(input$modtype, input$s.num, input$i.num,
        input$inf.prob, input$act.rate)
    isolate({
      if (is.null(rv$mod)) {
        withProgress(message = "Running initial model...", {
          rv$mod <- build_and_run()
        })
      }
    })
  })


  # =========================================================================
  # Summary
  # =========================================================================
  output$outSummary <- renderUI({
    m <- rv$mod
    req(m)

    p <- m$param
    nsims <- m$control$nsims
    type <- m$control$type
    nsteps <- m$control$nsteps

    # Get mean data for summary statistics
    df_mean <- as.data.frame(m, out = "mean")

    # Helper for styled stat boxes (matching DCM app)
    stat_box <- function(label, value, color = "#2C3E50") {
      div(
        class = "d-inline-block text-center me-4 mb-2",
        style = "min-width: 100px;",
        div(style = paste0("font-size: 1.17rem; font-weight: 500; color: ",
                           color, ";"),
            value),
        div(class = "text-muted", style = "font-size: 0.81rem;", label)
      )
    }

    # --- R0 / Force of Infection ---
    if (type == "SI") {
      foi <- p$inf.prob * p$act.rate
      r0_section <- tagList(
        tags$h6(class = "text-uppercase fw-semibold mb-2", "Transmission"),
        stat_box("Force of Infection", round(foi, 2), "#E74C3C"),
        p(class = "text-muted", style = "font-size: 0.81rem;",
          "In SI models, all infections are permanent.",
          "The force of infection (inf.prob \u00D7 act.rate) determines",
          "how quickly the susceptible population is depleted.")
      )
    } else {
      removal <- p$rec.rate
      if (!is.null(p$di.rate)) removal <- removal + p$di.rate
      r0 <- (p$inf.prob * p$act.rate) / removal
      if (r0 > 1) {
        r0_interp <- "greater than 1, so the infection is expected to spread."
      } else if (r0 == 1) {
        r0_interp <- "exactly 1, so the epidemic is at a tipping point."
      } else {
        r0_interp <- "less than 1, so the infection is expected to die out."
      }
      r0_section <- tagList(
        tags$h6(class = "text-uppercase fw-semibold mb-2", "Transmission"),
        stat_box(HTML("R<sub>0</sub>"), round(r0, 2),
                 if (r0 > 1) "#E74C3C" else "#18BC9C"),
        p(class = "text-muted", style = "font-size: 0.81rem;",
          HTML(paste0("R<sub>0</sub> is ", r0_interp,
                      " Each infected person generates on average <strong>",
                      round(r0, 2), "</strong> new infections",
                      " before recovering.")))
      )
    }

    # --- Peak & Timeline (from mean trajectory) ---
    peak_i <- max(df_mean$i.num)
    peak_t <- df_mean$time[which.max(df_mean$i.num)]
    peak_prev <- peak_i / df_mean$num[which.max(df_mean$i.num)]
    final_prev <- df_mean$i.num[nrow(df_mean)] / df_mean$num[nrow(df_mean)]

    # Cumulative infections and incidence metrics (from mean)
    has_vital <- !is.null(p$a.rate)
    attack_rate_valid <- type %in% c("SI", "SIR") && !has_vital
    if ("si.flow" %in% names(df_mean)) {
      cum_inf <- sum(df_mean$si.flow, na.rm = TRUE)
      cum_inc_rate <- cum_inf / sum(df_mean$s.num, na.rm = TRUE)
      if (attack_rate_valid) {
        attack_rate <- cum_inf / df_mean$s.num[1]
      }
    } else {
      cum_inf <- NA
      cum_inc_rate <- NA
    }

    timeline_section <- tagList(
      hr(),
      tags$h6(class = "text-uppercase fw-semibold mb-2",
              if (nsims > 1) "Epidemic Timeline (Mean)" else "Epidemic Timeline"),
      div(
        class = "d-flex flex-wrap",
        stat_box("Peak Infected", format(round(peak_i), big.mark = ","),
                 "#E74C3C"),
        stat_box("Peak Time", paste0("t = ", peak_t), "#3498DB"),
        stat_box("Peak Prevalence", paste0(round(peak_prev * 100, 1), "%"),
                 "#E74C3C")
      ),
      if (!is.na(cum_inf)) {
        div(
          class = "d-flex flex-wrap",
          stat_box("Cumulative Infections",
                   format(round(cum_inf), big.mark = ","), "#2C3E50"),
          stat_box("Incidence Rate",
                   paste0(round(cum_inc_rate * 1000, 2), " per 1k pt"),
                   "#2C3E50"),
          if (attack_rate_valid) {
            stat_box("Attack Rate",
                     paste0(round(attack_rate * 100, 1), "%"), "#2C3E50")
          }
        )
      },
      if (final_prev < 0.001 && type != "SI") {
        p(class = "text-muted", style = "font-size: 0.81rem;",
          "The epidemic has largely resolved by the end of the simulation.")
      } else if (type == "SI") {
        p(class = "text-muted", style = "font-size: 0.81rem;",
          "In SI models, prevalence increases monotonically toward 100%.")
      } else {
        p(class = "text-muted", style = "font-size: 0.81rem;",
          paste0("At the end of the simulation (t = ", nsteps,
                 "), mean prevalence is ", round(final_prev * 100, 1), "%."))
      }
    )

    # --- Stochastic Variation ---
    if (nsims > 1) {
      peaks <- numeric(nsims)
      peak_times <- numeric(nsims)
      final_prevs <- numeric(nsims)
      for (s in seq_len(nsims)) {
        df_s <- as.data.frame(m, out = "vals", sim = s)
        peaks[s] <- max(df_s$i.num, na.rm = TRUE)
        peak_times[s] <- df_s$time[which.max(df_s$i.num)]
        final_prevs[s] <- df_s$i.num[nrow(df_s)] / df_s$num[nrow(df_s)]
      }

      stoch_section <- tagList(
        hr(),
        tags$h6(class = "text-uppercase fw-semibold mb-2",
                "Stochastic Variation"),
        div(
          class = "d-flex flex-wrap",
          stat_box("Simulations", nsims, "#2C3E50"),
          stat_box("Peak Infected Range",
                   paste0(format(round(min(peaks)), big.mark = ","),
                          " \u2013 ",
                          format(round(max(peaks)), big.mark = ",")),
                   "#E74C3C"),
          stat_box("Peak Time Range",
                   paste0(min(peak_times), " \u2013 ", max(peak_times)),
                   "#3498DB")
        ),
        p(class = "text-muted", style = "font-size: 0.81rem;",
          paste0("Across ", nsims, " simulations, peak infected ranged from ",
                 format(round(min(peaks)), big.mark = ","), " to ",
                 format(round(max(peaks)), big.mark = ","),
                 " (SD = ", round(sd(peaks), 1), "). ",
                 "Final prevalence ranged from ",
                 round(min(final_prevs) * 100, 1), "% to ",
                 round(max(final_prevs) * 100, 1), "% ",
                 "(SD = ", round(sd(final_prevs) * 100, 1), " pp)."))
      )
    } else {
      stoch_section <- NULL
    }

    # --- Intervention impact ---
    if (!is.null(p$inter.eff) && !is.null(p$inter.start)) {
      inter_section <- tagList(
        hr(),
        tags$h6(class = "text-uppercase fw-semibold mb-2", "Intervention"),
        p(style = "font-size: 0.81rem;",
          paste0("An intervention reducing transmission by ",
                 round(p$inter.eff * 100), "% begins at t = ",
                 p$inter.start, ". ",
                 "Effective transmission probability drops from ",
                 round(p$inf.prob, 3), " to ",
                 round(p$inf.prob * (1 - p$inter.eff), 3),
                 " per act."))
      )
    } else {
      inter_section <- NULL
    }

    # --- Vital dynamics ---
    if (!is.null(p$a.rate) && p$vital) {
      # Build death rate description
      rates_same <- (p$ds.rate == p$di.rate) &&
        (is.null(p$dr.rate) || p$ds.rate == p$dr.rate)
      if (rates_same) {
        death_desc <- paste0("Deaths occur at rate ", p$ds.rate,
                             " across all compartments.")
      } else {
        death_parts <- paste0("S: ", p$ds.rate, ", I: ", p$di.rate)
        if (!is.null(p$dr.rate)) {
          death_parts <- paste0(death_parts, ", R: ", p$dr.rate)
        }
        death_desc <- paste0("Death rates per compartment: ", death_parts, ".")
      }
      vital_section <- tagList(
        hr(),
        tags$h6(class = "text-uppercase fw-semibold mb-2", "Vital Dynamics"),
        p(class = "text-muted", style = "font-size: 0.81rem;",
          paste0("Births enter at rate ", p$a.rate,
                 " per person per time step. ",
                 death_desc, " ",
                 "Final mean population size: ",
                 format(round(df_mean$num[nrow(df_mean)]), big.mark = ","),
                 " (started at ",
                 format(round(df_mean$num[1]), big.mark = ","), ")."))
      )
    } else {
      vital_section <- NULL
    }

    # --- Model configuration ---
    config_section <- tagList(
      hr(),
      tags$h6(class = "text-uppercase fw-semibold mb-2",
              "Model Configuration"),
      tags$table(
        class = "table table-sm table-borderless",
        style = "max-width: 400px; font-size: 0.81rem;",
        tags$tr(tags$td(class = "text-muted", "Model Class"),
                tags$td(tags$strong("ICM (Stochastic)"))),
        tags$tr(tags$td(class = "text-muted", "Disease Type"),
                tags$td(tags$strong(type))),
        tags$tr(tags$td(class = "text-muted", "Population"),
                tags$td(format(round(df_mean$num[1]), big.mark = ","))),
        tags$tr(tags$td(class = "text-muted", "Time Steps"),
                tags$td(nsteps)),
        tags$tr(tags$td(class = "text-muted", "Simulations"),
                tags$td(nsims)),
        tags$tr(tags$td(class = "text-muted", "Transmission Prob."),
                tags$td(p$inf.prob)),
        tags$tr(tags$td(class = "text-muted", "Act Rate"),
                tags$td(p$act.rate)),
        if (type != "SI") {
          tags$tr(tags$td(class = "text-muted", "Recovery Rate"),
                  tags$td(p$rec.rate))
        }
      )
    )

    # Assemble
    tagList(
      r0_section,
      timeline_section,
      stoch_section,
      inter_section,
      vital_section,
      config_section
    )
  })


  # =========================================================================
  # Main Plot (plotly — directly in UI, no renderUI wrapper)
  # =========================================================================
  output$MainPlotly <- plotly::renderPlotly({
    m <- rv$mod
    req(m)

    nsims <- m$control$nsims

    # Determine what to plot based on compsel
    if (input$compsel == "Compartment Prevalence") {
      y_cols <- c("s.num", "i.num")
      if (m$control$type == "SIR") y_cols <- c(y_cols, "r.num")
      is_prev <- TRUE
      y_label <- "Prevalence"
    } else if (input$compsel == "Compartment Size") {
      y_cols <- c("s.num", "i.num")
      if (m$control$type == "SIR") y_cols <- c(y_cols, "r.num")
      is_prev <- FALSE
      y_label <- "Number"
    } else if (input$compsel == "Disease Incidence") {
      y_cols <- "si.flow"
      is_prev <- FALSE
      y_label <- "New Infections"
    } else if (input$compsel == "Recovery Flow") {
      y_cols <- "ir.flow"
      is_prev <- FALSE
      y_label <- "Recoveries"
    } else if (input$compsel == "Re-susceptibility Flow") {
      y_cols <- "is.flow"
      is_prev <- FALSE
      y_label <- "Re-susceptible"
    } else if (input$compsel == "Arrival Flow") {
      y_cols <- "a.flow"
      is_prev <- FALSE
      y_label <- "Arrivals"
    } else if (input$compsel == "Departure Flows") {
      y_cols <- c("ds.flow", "di.flow")
      if (m$control$type == "SIR") y_cols <- c(y_cols, "dr.flow")
      is_prev <- FALSE
      y_label <- "Departures"
    } else {
      y_cols <- "i.num"
      is_prev <- FALSE
      y_label <- "Number"
    }

    # Nice labels and colors
    label_map <- c("s.num" = "Susceptible", "i.num" = "Infected",
                    "r.num" = "Recovered", "si.flow" = "S \u2192 I",
                    "ir.flow" = "I \u2192 R", "is.flow" = "I \u2192 S",
                    "a.flow" = "Arrivals",
                    "ds.flow" = "Deaths (S)", "di.flow" = "Deaths (I)",
                    "dr.flow" = "Deaths (R)")
    col_map <- c("s.num" = "#3498DB", "i.num" = "#E74C3C",
                  "r.num" = "#18BC9C", "si.flow" = "#E74C3C",
                  "ir.flow" = "#18BC9C", "is.flow" = "#3498DB",
                  "a.flow" = "#18BC9C",
                  "ds.flow" = "#3498DB", "di.flow" = "#E74C3C",
                  "dr.flow" = "#18BC9C")

    p <- plotly::plot_ly()

    if (nsims == 1) {
      # --- Single simulation: simple lines ---
      df <- as.data.frame(m, out = "vals")
      if (is_prev) {
        for (col in y_cols) {
          df[[paste0(col, ".prev")]] <- df[[col]] / df$num
        }
        plot_cols <- paste0(y_cols, ".prev")
      } else {
        plot_cols <- y_cols
      }
      for (i in seq_along(plot_cols)) {
        col <- plot_cols[i]
        orig <- y_cols[i]
        clr <- col_map[[orig]]
        lbl <- label_map[[orig]]
        p <- plotly::add_trace(p, x = df$time, y = df[[col]],
                               type = "scatter", mode = "lines",
                               name = lbl,
                               line = list(color = clr, width = 2.5))
      }
    } else {
      # --- Multiple simulations: mean lines + IQR ribbons ---
      df_mean <- as.data.frame(m, out = "mean")
      df_q25 <- as.data.frame(m, out = "qnt", qval = 0.25)
      df_q75 <- as.data.frame(m, out = "qnt", qval = 0.75)

      if (is_prev) {
        for (col in y_cols) {
          df_mean[[paste0(col, ".prev")]] <- df_mean[[col]] / df_mean$num
          df_q25[[paste0(col, ".prev")]] <- df_q25[[col]] / df_q25$num
          df_q75[[paste0(col, ".prev")]] <- df_q75[[col]] / df_q75$num
        }
        plot_cols <- paste0(y_cols, ".prev")
      } else {
        plot_cols <- y_cols
      }

      for (i in seq_along(plot_cols)) {
        col <- plot_cols[i]
        orig <- y_cols[i]
        clr <- col_map[[orig]]
        lbl <- label_map[[orig]]

        # IQR ribbon: lower bound (invisible line)
        p <- plotly::add_trace(p,
          x = df_q25$time, y = df_q25[[col]],
          type = "scatter", mode = "lines",
          line = list(color = "transparent"),
          showlegend = FALSE, name = paste0(lbl, " Q25"),
          hoverinfo = "skip")

        # IQR ribbon: upper bound (filled to lower)
        p <- plotly::add_trace(p,
          x = df_q75$time, y = df_q75[[col]],
          type = "scatter", mode = "lines",
          fill = "tonexty",
          fillcolor = paste0(clr, "33"),
          line = list(color = "transparent"),
          showlegend = FALSE, name = paste0(lbl, " IQR"),
          hoverinfo = "skip")

        # Mean line
        p <- plotly::add_trace(p,
          x = df_mean$time, y = df_mean[[col]],
          type = "scatter", mode = "lines",
          name = lbl,
          line = list(color = clr, width = 2.5))
      }
    }

    # Intervention line
    if (isTRUE(input$enable_intervention) &&
        !is.null(m$param$inter.eff)) {
      p <- plotly::layout(
        p,
        shapes = list(
          list(type = "line",
               x0 = m$param$inter.start, x1 = m$param$inter.start,
               y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "#95A5A6", width = 1.5, dash = "dash"))
        ),
        annotations = list(
          list(x = m$param$inter.start, y = 1, yref = "paper",
               text = "intervention", showarrow = FALSE,
               font = list(color = "#95A5A6", size = 11))
        )
      )
    }

    p <- plotly::layout(p,
                        xaxis = list(title = "Time"),
                        yaxis = list(title = y_label),
                        legend = list(x = 1.02, y = 1),
                        margin = list(t = 30))
    p <- plotly::config(p, displayModeBar = TRUE,
                        modeBarButtonsToRemove = c("lasso2d", "select2d"))
    p
  })


  # =========================================================================
  # Download main plot (base R plot.icm for PDF)
  # =========================================================================
  output$dlMainPlot <- downloadHandler(
    filename = "EpiModel_ICM_Plot.pdf",
    content = function(file) {
      m <- rv$mod
      req(m)
      pdf(file = file, height = 8, width = 12)
      par(mar = c(3.5, 3.5, 1.2, 1), mgp = c(2.1, 1, 0))
      if (input$compsel == "Compartment Prevalence") {
        plot(m, popfrac = TRUE, mean.line = TRUE, sim.lines = FALSE,
             qnts = 0.5, lwd = 3, legend = TRUE,
             leg.cex = 1.1, main = "")
      } else if (input$compsel == "Compartment Size") {
        plot(m, popfrac = FALSE, mean.line = TRUE, sim.lines = FALSE,
             qnts = 0.5, lwd = 3, legend = TRUE,
             leg.cex = 1.1, main = "")
      } else if (input$compsel == "Disease Incidence") {
        plot(m, y = "si.flow", popfrac = FALSE, mean.line = TRUE,
             sim.lines = FALSE, qnts = 0.5, lwd = 3,
             legend = TRUE, leg.cex = 1.1, main = "")
      } else {
        plot(m, popfrac = TRUE, mean.line = TRUE, sim.lines = FALSE,
             qnts = 0.5, lwd = 3, legend = TRUE,
             leg.cex = 1.1, main = "")
      }
      dev.off()
    }
  )


  # =========================================================================
  # Data Table
  # =========================================================================

  # Simulation number control for individual sim view
  output$simnoControl <- renderUI({
    m <- rv$mod
    req(m)
    maxsims <- m$control$nsims
    numericInput("datasim", "Simulation Number",
                 value = 1, min = 1, max = maxsims, step = 1)
  })

  output$outData <- DT::renderDT({
    m <- rv$mod
    req(m)
    if (input$datasel == "Means") {
      round(as.data.frame(m, out = "mean"), 2)
    } else if (input$datasel == "Standard Deviations") {
      round(as.data.frame(m, out = "sd"), 2)
    } else if (input$datasel == "Individual Simulations") {
      sim_num <- if (!is.null(input$datasim)) max(1, input$datasim) else 1
      as.data.frame(m, out = "vals", sim = sim_num)
    }
  }, options = list(lengthMenu = c(10, 25, 50, 100),
                    pageLength = 10,
                    scrollX = TRUE))

  output$dlData <- downloadHandler(
    filename = "EpiModel_ICM_Data.csv",
    content = function(file) {
      m <- rv$mod
      if (input$datasel == "Means") {
        write.csv(as.data.frame(m, out = "mean"), file, row.names = FALSE)
      } else if (input$datasel == "Standard Deviations") {
        write.csv(as.data.frame(m, out = "sd"), file, row.names = FALSE)
      } else if (input$datasel == "Individual Simulations") {
        sim_num <- if (!is.null(input$datasim)) max(1, input$datasim) else 1
        write.csv(as.data.frame(m, out = "vals", sim = sim_num),
                  file, row.names = FALSE)
      }
    }
  )

})
