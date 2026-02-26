##
## Server File for epidcm Shiny Application
##
## Run local: epiweb(class = "dcm")
##

library(shiny)
library(EpiModel)

# Compartment color palette
comp_colors <- c(
  "s.num" = "#3498DB",  # blue
  "i.num" = "#E74C3C",  # red
  "r.num" = "#18BC9C",  # green
  "num"   = "#95A5A6"   # gray
)

shinyServer(function(input, output, session) {

  # =========================================================================
  # Presets
  # =========================================================================
  observeEvent(input$preset, {
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
      updateSliderInput(session, "inf.prob", value = 0.9)
      updateSliderInput(session, "act.rate", value = 12)
      updateSliderInput(session, "rec.rate", value = round(1 / 8, 3))
      updateNumericInput(session, "s.num", value = 10000)
      updateNumericInput(session, "i.num", value = 1)
      updateNumericInput(session, "r.num", value = 0)
      updateNumericInput(session, "nsteps", value = 100)
      updateCheckboxInput(session, "enable_vital", value = FALSE)
      updateCheckboxInput(session, "enable_intervention", value = FALSE)
      updateCheckboxInput(session, "enable_sensitivity", value = FALSE)
    }
  })

  # =========================================================================
  # Dynamic UI updates
  # =========================================================================

  # Update flow choices based on model type and vital dynamics
  observe({
    type <- input$modtype
    vital <- input$enable_vital
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

  # Update time slider range when nsteps changes
  observe({
    ns <- input$nsteps
    if (!is.null(ns) && !is.na(ns) && ns >= 1) {
      updateSliderInput(session, "summTs",
                        max = ns,
                        value = min(input$summTs, ns))
    }
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
    if (input$enable_intervention) {
      p$inter.eff <- input$inter.eff
      p$inter.start <- input$inter.start
    }

    # Vital dynamics
    if (input$enable_vital) {
      p$a.rate <- input$a.rate
      p$ds.rate <- input$death_rate
      p$di.rate <- input$death_rate
      if (input$modtype == "SIR") {
        p$dr.rate <- input$death_rate
      }
    }

    # Sensitivity analysis: override selected param with vector
    if (input$enable_sensitivity) {
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
                odemethod = "rk4",
                verbose = FALSE)
  })

  mod <- eventReactive(input$runMod, {
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
  }, ignoreNULL = FALSE)


  # =========================================================================
  # Key Measures
  # =========================================================================
  output$keyMeasures <- renderUI({
    m <- mod()
    req(m)

    # R0 calculation
    p <- m$param
    if (input$modtype == "SI") {
      r0_val <- p$inf.prob[1] * p$act.rate[1]
      r0_label <- "Force of Infection"
      r0_formula <- "inf.prob \u00D7 act.rate"
    } else {
      removal <- p$rec.rate[1]
      if (!is.null(p$di.rate)) removal <- removal + p$di.rate[1]
      r0_val <- (p$inf.prob[1] * p$act.rate[1]) / removal
      r0_label <- HTML("R<sub>0</sub>")
      r0_formula <- "inf.prob \u00D7 act.rate / removal rate"
    }

    # Get time-specific measures (use run 1 for multi-run models)
    ts <- input$summTs
    df <- as.data.frame(m, run = 1)
    ts_row <- df[df$time == ts, ]
    if (nrow(ts_row) > 0) {
      row <- ts_row[1, ]
      prev <- round(row$i.num / row$num, 4)
      if ("si.flow" %in% names(row) && !is.na(row$si.flow)) {
        inc <- round(row$si.flow, 1)
      } else {
        inc <- NA
      }
    } else {
      prev <- NA
      inc <- NA
    }

    tagList(
      div(class = "mb-2",
          tags$strong(r0_label), ": ",
          tags$span(round(r0_val, 2), class = "text-primary fs-5"),
          br(),
          tags$small(class = "text-muted", r0_formula)
      ),
      hr(class = "my-2"),
      div(class = "mb-2",
          tags$strong("Prevalence"),
          tags$small(paste0(" (t=", ts, ")")), ": ",
          tags$span(
            if (!is.na(prev)) paste0(round(prev * 100, 1), "%") else "---",
            class = "text-danger"
          )
      ),
      div(class = "mb-1",
          tags$strong("New Infections"),
          tags$small(paste0(" (t=", ts, ")")), ": ",
          tags$span(
            if (!is.na(inc)) inc else "---",
            class = "text-warning"
          )
      )
    )
  })


  # =========================================================================
  # Main Plot (base R)
  # =========================================================================
  base_plot <- function() {
    m <- mod()
    req(m)

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
    } else if (input$compsel == "Recovery Flow") {
      plot(m, y = "ir.flow", popfrac = FALSE, lwd = 3,
           legend = "n", main = "")
    } else if (input$compsel == "Re-susceptibility Flow") {
      plot(m, y = "is.flow", popfrac = FALSE, lwd = 3,
           legend = "n", main = "")
    } else if (input$compsel == "Arrival Flow") {
      plot(m, y = "a.flow", popfrac = FALSE, lwd = 3,
           legend = "n", main = "")
    } else if (input$compsel == "Departure Flows") {
      y_flows <- "ds.flow"
      if ("di.flow" %in% names(m$epi)) y_flows <- c(y_flows, "di.flow")
      if ("dr.flow" %in% names(m$epi)) y_flows <- c(y_flows, "dr.flow")
      plot(m, y = y_flows, popfrac = FALSE, lwd = 3,
           legend = "full", leg.cex = 1.1, main = "")
    }

    # Draw intervention line
    if (input$enable_intervention && !is.null(m$param$inter.eff)) {
      abline(v = m$param$inter.start, lty = 2, col = "#95A5A6", lwd = 1.5)
      mtext("intervention", side = 3, at = m$param$inter.start,
            cex = 0.7, col = "#95A5A6")
    }
  }

  output$MainPlot <- renderPlot({
    base_plot()
  })


  # =========================================================================
  # Main Plot (plotly)
  # =========================================================================
  output$plotlyUI <- renderUI({
    if (requireNamespace("plotly", quietly = TRUE)) {
      plotly::plotlyOutput("MainPlotly", height = "450px")
    } else {
      helpText("Install the 'plotly' package for interactive plots: ",
               tags$code("install.packages(\"plotly\")"))
    }
  })

  if (requireNamespace("plotly", quietly = TRUE)) {
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

      # Color map
      col_map <- c("s.num" = "#3498DB", "s.num.prev" = "#3498DB",
                    "i.num" = "#E74C3C", "i.num.prev" = "#E74C3C",
                    "r.num" = "#18BC9C", "r.num.prev" = "#18BC9C",
                    "si.flow" = "#E74C3C", "ir.flow" = "#18BC9C",
                    "is.flow" = "#3498DB",
                    "a.flow" = "#18BC9C",
                    "ds.flow" = "#3498DB", "di.flow" = "#E74C3C",
                    "dr.flow" = "#18BC9C")

      # Build plot
      if (nruns == 1) {
        p <- plotly::plot_ly()
        for (col in plot_cols) {
          clr <- ifelse(col %in% names(col_map), col_map[[col]], "#2C3E50")
          p <- plotly::add_trace(p, x = df$time, y = df[[col]],
                                 type = "scatter", mode = "lines",
                                 name = col,
                                 line = list(color = clr, width = 2.5))
        }
      } else {
        p <- plotly::plot_ly()
        run_names <- names(m$epi[[plot_cols[1]]])
        palette <- grDevices::colorRampPalette(c("#3498DB", "#E74C3C"))(nruns)
        for (i in seq_len(nruns)) {
          for (col in plot_cols) {
            run_df <- as.data.frame(m, run = i)
            p <- plotly::add_trace(
              p, x = run_df$time, y = run_df[[col]],
              type = "scatter", mode = "lines",
              name = if (length(plot_cols) == 1) {
                run_names[i]
              } else {
                paste(col, run_names[i])
              },
              line = list(color = palette[i], width = 2),
              legendgroup = run_names[i]
            )
          }
        }
      }

      # Intervention line
      if (input$enable_intervention && !is.null(m$param$inter.eff)) {
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
  }


  # =========================================================================
  # Download main plot
  # =========================================================================
  output$dlMainPlot <- downloadHandler(
    filename = "EpiModel_DCM_Plot.pdf",
    content = function(file) {
      pdf(file = file, height = 6, width = 10)
      base_plot()
      dev.off()
    }
  )


  # =========================================================================
  # Compartment Diagram
  # =========================================================================
  output$CompPlot <- renderPlot({
    m <- mod()
    req(m)

    ts <- input$summTs
    if (is.na(ts)) ts <- 1
    ts <- max(1, min(ts, m$control$nsteps))

    # comp_plot only supports 1-group models
    if (m$param$groups == 1) {
      comp_plot(m, at = ts, digits = 1)
    }
  })


  # =========================================================================
  # Summary
  # =========================================================================
  output$outSummary <- renderPrint({
    m <- mod()
    req(m)

    ts <- input$summTs
    if (is.na(ts)) ts <- 1
    ts <- max(1, min(ts, m$control$nsteps))

    summary(m, at = ts, digits = 3)
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
