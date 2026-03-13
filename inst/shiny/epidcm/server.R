##
## Server File for epidcm Shiny Application
##
## Run local: epiweb(class = "dcm")
##

library(shiny)
library(EpiModel)

shinyServer(function(input, output, session) {

  # =========================================================================
  # Presets
  # =========================================================================
  preset_updating <- reactiveVal(FALSE)

  observeEvent(input$preset, {
    preset_updating(TRUE)
    on.exit(preset_updating(FALSE))
    if (input$preset == "Flu-like (SIR)") {
      updateSelectInput(session, "modtype", selected = "SIR")
      updateSliderInput(session, "inf.prob", value = 0.03)
      updateSliderInput(session, "act.rate", value = 10)
      updateSliderInput(session, "rec.rate", value = round(1 / 7, 3))
      updateNumericInput(session, "s.num", value = 10000)
      updateNumericInput(session, "i.num", value = 1)
      updateNumericInput(session, "r.num", value = 0)
      updateNumericInput(session, "nsteps", value = 200)
      updateCheckboxInput(session, "enable_vital", value = FALSE)
      updateCheckboxInput(session, "enable_intervention", value = FALSE)
      updateCheckboxInput(session, "enable_sensitivity", value = FALSE)
    } else if (input$preset == "STI-like (SIS)") {
      updateSelectInput(session, "modtype", selected = "SIS")
      updateSliderInput(session, "inf.prob", value = 0.2)
      updateSliderInput(session, "act.rate", value = 0.5)
      updateSliderInput(session, "rec.rate", value = 0.01)
      updateNumericInput(session, "s.num", value = 1000)
      updateNumericInput(session, "i.num", value = 10)
      updateNumericInput(session, "nsteps", value = 500)
      updateCheckboxInput(session, "enable_vital", value = FALSE)
      updateCheckboxInput(session, "enable_intervention", value = FALSE)
      updateCheckboxInput(session, "enable_sensitivity", value = FALSE)
    } else if (input$preset == "Measles-like (SIR)") {
      updateSelectInput(session, "modtype", selected = "SIR")
      updateSliderInput(session, "inf.prob", value = 0.5)
      updateSliderInput(session, "act.rate", value = 3)
      updateSliderInput(session, "rec.rate", value = 0.1)
      updateNumericInput(session, "s.num", value = 10000)
      updateNumericInput(session, "i.num", value = 1)
      updateNumericInput(session, "r.num", value = 0)
      updateNumericInput(session, "nsteps", value = 100)
      updateCheckboxInput(session, "enable_vital", value = FALSE)
      updateCheckboxInput(session, "enable_intervention", value = FALSE)
      updateCheckboxInput(session, "enable_sensitivity", value = FALSE)
    }
  })

  # Reset preset to "Custom" when user manually changes disease type
  observeEvent(input$modtype, {
    if (!preset_updating() && input$preset != "Custom") {
      updateSelectInput(session, "preset", selected = "Custom")
    }
  })

  # =========================================================================
  # Dynamic UI updates
  # =========================================================================

  # Update flow choices based on model type and vital dynamics
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

  # Update sensitivity parameter choices based on model type
  observe({
    type <- input$modtype
    choices <- c("inf.prob", "act.rate")
    if (type != "SI") {
      choices <- c(choices, "rec.rate")
    }
    updateSelectInput(session, "sens_param", choices = choices)
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
  # Core model reactive
  # =========================================================================
  param <- reactive({
    p <- list(inf.prob = input$inf.prob,
              act.rate = input$act.rate)

    # Recovery rate
    if (input$modtype != "SI") {
      p$rec.rate <- input$rec.rate
    }

    # Intervention
    if (isTRUE(input$enable_intervention)) {
      p$inter.eff <- input$inter.eff
      p$inter.start <- input$inter.start
    }

    # Vital dynamics
    if (isTRUE(input$enable_vital)) {
      p$a.rate <- input$a.rate
      p$ds.rate <- input$ds.rate
      if (isTRUE(input$diff_death_rates)) {
        p$di.rate <- input$di.rate
        if (input$modtype == "SIR") {
          p$dr.rate <- input$dr.rate
        }
      } else {
        p$di.rate <- input$ds.rate
        if (input$modtype == "SIR") {
          p$dr.rate <- input$ds.rate
        }
      }
    }

    # Sensitivity analysis: override selected param with vector
    if (isTRUE(input$enable_sensitivity)) {
      nruns <- as.integer(input$sens_nruns)
      sens_vals <- seq(input$sens_min, input$sens_max, length.out = nruns)
      p[[input$sens_param]] <- sens_vals
    }

    do.call(param.dcm, p)
  })

  init <- reactive({
    args <- list(s.num = input$s.num, i.num = input$i.num)
    if (input$modtype == "SIR") {
      args$r.num <- input$r.num
    }
    do.call(init.dcm, args)
  })

  control <- reactive({
    control.dcm(type = input$modtype,
                nsteps = input$nsteps,
                dt = 1,
                odemethod = "lsoda",
                verbose = FALSE)
  })

  mod <- reactive({
    validate(
      need(input$s.num > 0, "Number susceptible must be positive"),
      need(input$i.num >= 0, "Number infected cannot be negative"),
      need(input$nsteps >= 10, "Need at least 10 time steps"),
      need(input$inf.prob > 0 && input$inf.prob <= 1,
           "Transmission probability must be between 0 and 1"),
      need(input$act.rate >= 0, "Act rate cannot be negative")
    )
    if (input$modtype != "SI") {
      validate(need(input$rec.rate >= 0, "Recovery rate cannot be negative"))
    }
    tryCatch(
      dcm(param(), init(), control()),
      error = function(e) {
        showNotification(paste("Model error:", e$message),
                         type = "error", duration = 8)
        NULL
      }
    )
  })


  # =========================================================================
  # Main Plot (plotly — directly in UI, no renderUI wrapper)
  # =========================================================================
  output$MainPlotly <- plotly::renderPlotly({
    m <- mod()
    req(m)

    df <- as.data.frame(m)
    nruns <- m$control$nruns

    # Determine what to plot
    if (input$compsel == "Compartment Prevalence") {
      y_cols <- c("s.num", "i.num")
      if (input$modtype == "SIR") y_cols <- c(y_cols, "r.num")
      for (col in y_cols) {
        df[[paste0(col, ".prev")]] <- df[[col]] / df$num
      }
      plot_cols <- paste0(y_cols, ".prev")
      y_label <- "Prevalence"
    } else if (input$compsel == "Compartment Size") {
      plot_cols <- c("s.num", "i.num")
      if (input$modtype == "SIR") plot_cols <- c(plot_cols, "r.num")
      y_label <- "Number"
    } else if (input$compsel == "Disease Incidence") {
      plot_cols <- "si.flow"
      y_label <- "New Infections"
    } else if (input$compsel == "Recovery Flow") {
      plot_cols <- "ir.flow"
      y_label <- "Recoveries"
    } else if (input$compsel == "Re-susceptibility Flow") {
      plot_cols <- "is.flow"
      y_label <- "Re-susceptible"
    } else if (input$compsel == "Arrival Flow") {
      plot_cols <- "a.flow"
      y_label <- "Arrivals"
    } else if (input$compsel == "Departure Flows") {
      plot_cols <- intersect(c("ds.flow", "di.flow", "dr.flow"), names(df))
      y_label <- "Departures"
    } else {
      plot_cols <- "i.num"
      y_label <- "Number"
    }

    # Color and label maps
    col_map <- c("s.num" = "#3498DB", "s.num.prev" = "#3498DB",
                  "i.num" = "#E74C3C", "i.num.prev" = "#E74C3C",
                  "r.num" = "#18BC9C", "r.num.prev" = "#18BC9C",
                  "si.flow" = "#E74C3C", "ir.flow" = "#18BC9C",
                  "is.flow" = "#3498DB",
                  "a.flow" = "#18BC9C",
                  "ds.flow" = "#3498DB", "di.flow" = "#E74C3C",
                  "dr.flow" = "#18BC9C")
    label_map <- c("s.num" = "Susceptible", "s.num.prev" = "Susceptible",
                    "i.num" = "Infected", "i.num.prev" = "Infected",
                    "r.num" = "Recovered", "r.num.prev" = "Recovered",
                    "si.flow" = "S \u2192 I",
                    "ir.flow" = "I \u2192 R", "is.flow" = "I \u2192 S",
                    "a.flow" = "Arrivals",
                    "ds.flow" = "Deaths (S)", "di.flow" = "Deaths (I)",
                    "dr.flow" = "Deaths (R)")

    # Build plot
    if (nruns == 1) {
      p <- plotly::plot_ly()
      for (col in plot_cols) {
        clr <- ifelse(col %in% names(col_map), col_map[[col]], "#2C3E50")
        lbl <- ifelse(col %in% names(label_map), label_map[[col]], col)
        p <- plotly::add_trace(p, x = df$time, y = df[[col]],
                               type = "scatter", mode = "lines",
                               name = lbl,
                               line = list(color = clr, width = 2.5))
      }
    } else {
        p <- plotly::plot_ly()

        # Get sensitivity parameter values for legend labels
        sens_param_name <- input$sens_param
        sens_vals <- m$param[[sens_param_name]]
        run_labels <- paste0(sens_param_name, " = ", round(sens_vals, 4))

        # For multi-run, reduce to primary variable to match base R plot.dcm
        is_prev <- input$compsel == "Compartment Prevalence"
        if (is_prev) {
          plot_cols <- "i.num.prev"
          y_label <- "Prevalence (Infected)"
        } else if (input$compsel == "Compartment Size") {
          plot_cols <- "i.num"
          y_label <- "Number Infected"
        }

        palette <- grDevices::colorRampPalette(c("#3498DB", "#E74C3C"))(nruns)
        for (i in seq_len(nruns)) {
          run_df <- as.data.frame(m, run = i)
          if (is_prev) {
            run_df[["i.num.prev"]] <- run_df$i.num / run_df$num
          }
          for (col in plot_cols) {
            p <- plotly::add_trace(
              p, x = run_df$time, y = run_df[[col]],
              type = "scatter", mode = "lines",
              name = run_labels[i],
              line = list(color = palette[i], width = 2)
            )
          }
        }
      }

      # Intervention line
      if (isTRUE(input$enable_intervention) && !is.null(m$param$inter.eff)) {
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
  # Download main plot
  # =========================================================================
  output$dlMainPlot <- downloadHandler(
    filename = "EpiModel_DCM_Plot.pdf",
    content = function(file) {
      m <- mod()
      req(m)
      pdf(file = file, height = 8, width = 12)
      par(mar = c(3.5, 3.5, 1.2, 1), mgp = c(2.1, 1, 0))
      if (input$compsel == "Compartment Prevalence") {
        plot(m, popfrac = TRUE, lwd = 3, legend = "full",
             leg.cex = 1.1, main = "")
      } else if (input$compsel == "Compartment Size") {
        plot(m, popfrac = FALSE, lwd = 3, legend = "full",
             leg.cex = 1.1, main = "")
      } else if (input$compsel == "Disease Incidence") {
        plot(m, y = "si.flow", popfrac = FALSE, lwd = 3,
             legend = "n", main = "")
      } else {
        plot(m, popfrac = TRUE, lwd = 3, legend = "full",
             leg.cex = 1.1, main = "")
      }
      dev.off()
    }
  )


  # =========================================================================
  # Summary
  # =========================================================================
  output$outSummary <- renderUI({
    m <- mod()
    req(m)

    p <- m$param
    nruns <- m$control$nruns
    type <- m$control$type
    df <- as.data.frame(m, run = 1)
    nsteps <- m$control$nsteps

    # Helper for styled stat boxes
    stat_box <- function(label, value, color = "#2C3E50") {
      div(
        class = "d-inline-block text-center me-4 mb-2",
        style = "min-width: 100px;",
        div(style = paste0("font-size: 1.17rem; font-weight: 500; color: ", color, ";"),
            value),
        div(class = "text-muted", style = "font-size: 0.81rem;", label)
      )
    }

    # --- R0 / Force of Infection ---
    if (type == "SI") {
      foi <- p$inf.prob[1] * p$act.rate[1]
      r0_section <- tagList(
        tags$h6(class = "text-uppercase fw-semibold mb-2", "Transmission"),
        stat_box("Force of Infection", round(foi, 2), "#E74C3C"),
        p(class = "text-muted", style = "font-size: 0.81rem;",
          "In SI models, all infections are permanent.",
          "The force of infection (inf.prob \u00D7 act.rate) determines",
          "how quickly the susceptible population is depleted.")
      )
    } else {
      removal <- p$rec.rate[1]
      if (!is.null(p$di.rate)) removal <- removal + p$di.rate[1]
      r0 <- (p$inf.prob[1] * p$act.rate[1]) / removal
      if (r0 > 1) {
        r0_interp <- "greater than 1, so the infection will spread in the population."
      } else if (r0 == 1) {
        r0_interp <- "exactly 1, so the epidemic is at a tipping point."
      } else {
        r0_interp <- "less than 1, so the infection will die out."
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

    # --- Peak & Timeline (run 1) ---
    peak_i <- max(df$i.num)
    peak_t <- df$time[which.max(df$i.num)]
    peak_prev <- peak_i / df$num[which.max(df$i.num)]
    final_prev <- df$i.num[nrow(df)] / df$num[nrow(df)]

    # Cumulative infections
    if ("si.flow" %in% names(df)) {
      cum_inf <- sum(df$si.flow, na.rm = TRUE)
      attack_rate <- cum_inf / df$s.num[1]
    } else {
      cum_inf <- NA
      attack_rate <- NA
    }

    timeline_section <- tagList(
      hr(),
      tags$h6(class = "text-uppercase fw-semibold mb-2", "Epidemic Timeline"),
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
          stat_box("Attack Rate", paste0(round(min(attack_rate, 1) * 100, 1), "%"),
                   "#2C3E50")
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
                 "), prevalence is ", round(final_prev * 100, 1), "%."))
      }
    )

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
                 round(p$inf.prob[1], 3), " to ",
                 round(p$inf.prob[1] * (1 - p$inter.eff), 3),
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
                 "Final population size: ",
                 format(round(df$num[nrow(df)]), big.mark = ","),
                 " (started at ",
                 format(round(df$num[1]), big.mark = ","), ")."))
      )
    } else {
      vital_section <- NULL
    }

    # --- Sensitivity note ---
    if (nruns > 1) {
      sens_section <- tagList(
        hr(),
        tags$h6(class = "text-uppercase fw-semibold mb-2", "Sensitivity Analysis"),
        p(class = "text-muted", style = "font-size: 0.81rem;",
          paste0(nruns, " model runs varying ", input$sens_param,
                 " from ", min(p[[input$sens_param]]),
                 " to ", max(p[[input$sens_param]]), ".",
                 " Statistics above are shown for run 1."))
      )
    } else {
      sens_section <- NULL
    }

    # --- Model configuration ---
    config_section <- tagList(
      hr(),
      tags$h6(class = "text-uppercase fw-semibold mb-2", "Model Configuration"),
      tags$table(
        class = "table table-sm table-borderless",
        style = "max-width: 400px; font-size: 0.81rem;",
        tags$tr(tags$td(class = "text-muted", "Model Class"),
                tags$td(tags$strong("DCM (Deterministic)"))),
        tags$tr(tags$td(class = "text-muted", "Disease Type"),
                tags$td(tags$strong(type))),
        tags$tr(tags$td(class = "text-muted", "Population"),
                tags$td(format(round(df$num[1]), big.mark = ","))),
        tags$tr(tags$td(class = "text-muted", "Time Steps"),
                tags$td(nsteps)),
        tags$tr(tags$td(class = "text-muted", "Transmission Prob."),
                tags$td(p$inf.prob[1])),
        tags$tr(tags$td(class = "text-muted", "Act Rate"),
                tags$td(p$act.rate[1])),
        if (type != "SI") {
          tags$tr(tags$td(class = "text-muted", "Recovery Rate"),
                  tags$td(p$rec.rate[1]))
        }
      )
    )

    # Assemble
    tagList(
      r0_section,
      timeline_section,
      inter_section,
      vital_section,
      sens_section,
      config_section
    )
  })


  # =========================================================================
  # Data Table
  # =========================================================================
  output$outData <- DT::renderDT({
    m <- mod()
    req(m)
    round(as.data.frame(m), 2)
  }, options = list(lengthMenu = c(10, 25, 50, 100),
                    pageLength = 10,
                    scrollX = TRUE))

  output$dlData <- downloadHandler(
    filename = "EpiModel_DCM_Data.csv",
    content = function(file) {
      write.csv(as.data.frame(mod()), file, row.names = FALSE)
    }
  )

})
