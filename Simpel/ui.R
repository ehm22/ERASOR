library(DT)
library(shiny)
library(shinyjs)
library(shinythemes)
##&*## purple design
library(bslib)
##&*## purple design

# if this is TRUE the patient specific tab will show, if FALSE then the tab is hidden
show_patient_tab <- TRUE 
# dev mode keeps interface clickabel and usable so ui changes can be checked should be false for users since they will just use it and then get disappointed it doesnt work 
dev_mode <- FALSE

### for sliders with numeric input boxes
rangeFilterUI <- function(
    id,
    label = NULL,
    min,
    max,
    value,
    step = 1,
    fixed = c("none", "left", "right"),
    label_from = "From:",
    label_to = "To:",
    label_left = "Keep ASOs with a value lower than:",
    label_right = "Greater or equal to:"
) {
  
  fixed <- match.arg(fixed)
  ns <- NS(id)
  
  tagList(
    if (!is.null(label)) h4(label),
    
    if (fixed == "none") {
      tagList(
        fluidRow(
          column(
            6,
            numericInput(ns("min"), label_from, value = value[1], min = min, max = max, step = step)
          ),
          column(
            6,
            numericInput(ns("max"), label_to, value = value[2], min = min, max = max, step = step)
          )
        ),
        sliderInput(ns("slider"), label = NULL, min = min, max = max, value = value, step = step)
      )
      
    } else if (fixed == "left") {
      tagList(
        numericInput(ns("max"), label_left, value = value, min = min, max = max, step = step),
        sliderInput(ns("slider"), label = NULL, min = min, max = max, value = value, step = step)
      )
      
    } else { # fixed == "right"
      tagList(
        numericInput(ns("min"), label_right, value = value, min = min, max = max, step = step),
        div(
          class = "slider-right-fixed",
          sliderInput(ns("slider"), label = NULL, min = min, max = max, value = value, step = step)
        )
      )
    }
  )
}


filterBlockUI <- function(title, checkbox_id, checkbox_value = TRUE, tooltip = NULL, content) {
  tagList(
    div(
      class = "filter-block",
      
      div(
        class = "filter-title-row",
        h5(
          tagList(
            HTML(paste0("<b>", title, "</b>")),
            if (!is.null(tooltip)) {
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = tooltip,
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            }
          )
        ),
        checkboxInput(checkbox_id, "Enable", value = checkbox_value)
      ),
      
      div(
        class = "filter-control-area",
        content
      )
    )
  )
}



ui <- fluidPage(
  useShinyjs(),
  
  if (!dev_mode) {
    div(
      id = "app_startup_loading",
      div(
        id = "app_startup_loading_box",
        span(id = "startup_loading_text", "Please wait. Loading"),
        span(id = "loading_dots", ".")
      )
    )
  },
  
     theme = shinythemes::shinytheme("flatly"),
  
  #### purple design ##&*####
  # theme = bs_theme(
  #   bootswatch = "flatly",
  #   primary = "purple",
  #   bg = "#ffe4ec",   # light pink background
  #   fg = "#2b002b"    # dark text so it stays readable
  # ),
  #### purple design ##&*####
  
  tags$head(
    tags$style(
      HTML("
        #app_startup_loading {
          position: fixed;
          inset: 0;
          background: rgba(255, 255, 255, 0.35);
          z-index: 99999;
          display: flex;
          align-items: center;
          justify-content: center;
        }
        
        table.dataTable thead th {
          white-space: nowrap !important;
          vertical-align: middle !important;
          text-align: center !important;
        }
        
        table.dataTable thead th .col-header-with-tooltip {
          display: inline-flex;
          align-items: center;
          gap: 4px;
          white-space: nowrap;
        }
        
       .col-header-with-tooltip {
          display: inline-flex;
          align-items: center;
          gap: 4px;
          position: relative;
       }
        
        #aso_seq_input {
          resize: none;
        }
                
        .tooltip-icon {
          position: relative;
          display: inline-flex;
          align-items: center;
          cursor: help;
        }
        
        /* Main top-level tab titles */
        #mode_tabs_wrapper > .tabbable > ul.nav-tabs > li > a {
          font-size: 18px !important;
          font-weight: 600 !important;
        }   
        
        
        .custom-tooltip {
          visibility: hidden;
          opacity: 0;
          position: absolute;
          bottom: 125%;          /* place above icon */
          left: 50%;
          transform: translateX(-50%);
          background: #222;
          color: #fff;
          padding: 6px 8px;
          border-radius: 4px;
          font-size: 12px;
          max-width: 400px;
          width: 400px;
          white-space: normal;
          text-align: center;
          z-index: 10000;
          transition: opacity 0.15s ease;
        }
        
        /* show tooltip */
        .tooltip-icon:hover .custom-tooltip {
          visibility: visible;
          opacity: 1;
        }
                
        /* rows clickable */
          table.dataTable tbody tr {
          cursor: pointer;
          transition: background-color 0.15s ease;
        }

        /* hover highlight */
          table.dataTable tbody tr:hover {
          background-color: #e8f4ff !important;
        }
          
        /* file upload progress bar */
          .file-upload-complete .progress-bar {
            background-color: #5cb85c !important;
            border-color: #4cae4c !important;
        }
                  
        /* center DT column headers */
            table.dataTable thead th {
           text-align: center !important;
            }
            
        #toggle_cols:disabled,
        #toggle_cols.disabled,
        #toggle_cols_patient:disabled,
        #toggle_cols_patient.disabled,
        #toggle_cols_patient_snv:disabled,
        #toggle_cols_patient_snv.disabled,
        #toggle_cols_single_aso:disabled,
        #toggle_cols_single_aso.disabled,
        #Download_single_aso_results:disabled,
        #Download_single_aso_results.disabled,
        #single_aso_download_rnaseh:disabled,
        #single_aso_download_rnaseh.disabled,
        #single_aso_download_offtarget:disabled,
        #single_aso_download_offtarget.disabled {
          opacity: 0.65 !important;
          cursor: default !important;
        }
        
        #Download_single_aso_results.disabled,
        #single_aso_download_rnaseh.disabled,
        #single_aso_download_offtarget.disabled,
        #Download_single_aso_results[disabled],
        #single_aso_download_rnaseh[disabled],
        #single_aso_download_offtarget[disabled] {
          cursor: default !important;
          pointer-events: auto !important;
        }
        
        /* Make main Run, Set to default, and file Browse buttons darker */
        #run_button,
        #reset_defaults,
        #run_button_patient,
        #reset_defaults_patient,
        #run_button_single_aso,
        .btn-file {
          background-color: #2c3e50 !important;
          border-color: #1a252f !important;
          color: #ffffff !important;
        }
        
        /* Slightly wider sidebar than Shiny width = 3 */
        .wide-sidebar-layout .col-sm-3 {
          width: 27.5% !important;
        }
        
        .wide-sidebar-layout .col-sm-9 {
          width: 72.5% !important;
        }
        
                
        /* Space inside grey sidebar */
        .wide-sidebar-layout .well {
          padding-left: 20px !important;
          padding-right: 20px !important;
        }
        
        /* Shared horizontal alignment for ASO length, filter titles, numeric boxes, and sliders */
        .wide-sidebar-layout #oligo_length_block,
        .wide-sidebar-layout .filter-block {
          padding-left: 8px !important;
          padding-right: 8px !important;
          box-sizing: border-box !important;
        }
        
        /* Patient ASO length slider too */
        .wide-sidebar-layout .shiny-input-container:has(#patient_oligo_length_range) {
          padding-left: 8px !important;
          padding-right: 8px !important;
          box-sizing: border-box !important;
        }
        
        /* More space between filters */
        .wide-sidebar-layout .filter-block {
          margin-top: 30px !important;
          margin-bottom: 30px !important;
        }
        
        /* Filter title left, Enable checkbox right */
        .wide-sidebar-layout .filter-title-row {
          display: flex !important;
          align-items: center !important;
          justify-content: space-between !important;
          gap: 12px !important;
          margin-bottom: 6px !important;
        }
        
        /* Filter title styling */
        .wide-sidebar-layout .filter-title-row h5 {
          margin: 0 !important;
          font-size: 17px !important;
          line-height: 1.3 !important;
        }
        
        /* Enable checkbox on the right, same height as title */
        .wide-sidebar-layout .filter-title-row .checkbox,
        .wide-sidebar-layout .filter-title-row .form-group {
          margin: 0 !important;
        }
        
        /* Enable checkbox text*/
        .wide-sidebar-layout .filter-title-row label {
          margin-bottom: 0 !important;
          font-size: 15px !important;
          font-weight: normal !important;
          white-space: nowrap !important;
        }
        
        /* Important: no extra nested side padding, so numeric inputs and sliders align with filter title */
        .wide-sidebar-layout .filter-control-area {
          padding-left: 0 !important;
          padding-right: 0 !important;
          box-sizing: border-box !important;
        }
        
        /* Bootstrap rows normally have negative margins; remove them inside filter blocks */
        .wide-sidebar-layout .filter-control-area .row,
        .wide-sidebar-layout .filter-control-area .fluid-row {
          margin-left: 0 !important;
          margin-right: 0 !important;
        }
        
        /* Numeric input columns should not create extra outer offset */
        .wide-sidebar-layout .filter-control-area .col-sm-6,
        .wide-sidebar-layout .filter-control-area .col-md-6,
        .wide-sidebar-layout .filter-control-area .col-lg-6,
        .wide-sidebar-layout .filter-control-area .col-xs-6 {
          padding-left: 0 !important;
          padding-right: 8px !important;
        }
        
        .wide-sidebar-layout .filter-control-area .col-sm-6:last-child,
        .wide-sidebar-layout .filter-control-area .col-md-6:last-child,
        .wide-sidebar-layout .filter-control-area .col-lg-6:last-child,
        .wide-sidebar-layout .filter-control-area .col-xs-6:last-child {
          padding-left: 8px !important;
          padding-right: 0 !important;
        }
        
        /* Numeric input labels/
        .wide-sidebar-layout .filter-control-area .form-group > label.control-label {
          font-size: 15px !important;
          font-weight: 400 !important;
          margin-left: 0 !important;
          padding-left: 0 !important;
        }
        
        /* Numeric input labels: normal weight and slightly smaller */
        .wide-sidebar-layout .filter-control-area .shiny-input-container label,
        .wide-sidebar-layout .filter-control-area .shiny-input-container label.control-label,
        .wide-sidebar-layout .filter-control-area .shiny-input-container label span,
        .wide-sidebar-layout .filter-control-area .shiny-input-container .control-label {
          font-size: 14px !important;
          font-weight: 400 !important;
          font-family: inherit !important;
        }
        
        /* Numeric input boxes aligned with filter title */
        .wide-sidebar-layout .filter-control-area input[type='number'] {
          width: 100% !important;
          box-sizing: border-box !important;
        }
        
        /* Match ASO length title size and weight with filter titles */
        .wide-sidebar-layout #oligo_length_block .control-label,
        .wide-sidebar-layout #patient_oligo_length_block .control-label {
          font-size: 17px !important;
          line-height: 1.3 !important;
          font-weight: 700 !important;
        }
        
        /* Slider aligned with filter title and numeric boxes */
        .wide-sidebar-layout .filter-control-area .shiny-input-container:has(.irs) {
          width: 100% !important;
          margin-top: 8px !important;
          margin-bottom: 6px !important;
          padding-left: 0 !important;
          padding-right: 0 !important;
          box-sizing: border-box !important;
        }
        
        /* Slider size */
        .wide-sidebar-layout .irs {
          height: 70px !important;
          margin-top: 6px !important;
          margin-bottom: 6px !important;
        }
        
        /* If the slider has grid/tick labels, give it room */
        .wide-sidebar-layout .irs-with-grid {
          height: 84px !important;
        }
        
        /* Hide min and max labels above sliders */
        .wide-sidebar-layout .irs-min,
        .wide-sidebar-layout .irs-max {
          display: none !important;
        }
        
        /* Slightly smaller current value bubble */
        .wide-sidebar-layout .irs-single,
        .wide-sidebar-layout .irs-from,
        .wide-sidebar-layout .irs-to {
          top: 0px !important;
          font-size: 14px !important;
          padding: 2px 6px !important;
          line-height: 1.25 !important;
          min-width: 26px !important;
          text-align: center !important;
          white-space: nowrap !important;
        }
        
        /* Track line */
        .wide-sidebar-layout .irs-line {
          top: 34px !important;
        }
        
        /* Filled bar */
        .wide-sidebar-layout .irs-bar {
          top: 34px !important;
        }
        
        /* Slider handle/ball */
        .wide-sidebar-layout .irs-handle {
          top: 26px !important;
        }
        
        /* Tick/grid labels */
        .wide-sidebar-layout .irs-grid {
          top: 52px !important;
        }
        
        /* Slightly smaller tick labels */
        .wide-sidebar-layout .irs-grid-text {
          font-size: 13.5px !important;
        }
        
        /* Keep minor ticks hidden where requested */
        .no-minor-ticks .irs-grid-pol.small {
          display: none !important;
        }

        /* Hover state */
        #run_button:hover,
        #reset_defaults:hover,
        #run_button_patient:hover,
        #reset_defaults_patient:hover,
        #run_button_single_aso:hover,
        .btn-file:hover {
          background-color: #1a252f !important;
          border-color: #101820 !important;
          color: #ffffff !important;
        }
        
        /* Active/clicked state */
        #run_button:active,
        #reset_defaults:active,
        #run_button_patient:active,
        #reset_defaults_patient:active,
        #run_button_single_aso:active,
        .btn-file:active {
          background-color: #101820 !important;
          border-color: #0b1117 !important;
          color: #ffffff !important;
        }
        
        .disabled-section,
        .disabled-section * {
          color: #999999 !important;
        }
        
        .disabled-section {
          opacity: 0.55 !important;
          pointer-events: none;
        }
        
        .disabled-section .irs-bar,
        .disabled-section .irs-line,
        .disabled-section .irs-handle,
        .disabled-section .irs-single,
        .disabled-section .irs-from,
        .disabled-section .irs-to {
          opacity: 0.55 !important;
        }
        
        .disabled-section .control-label,
        .disabled-section h4,
        .disabled-section h5,
        .disabled-section label,
        .disabled-section .shiny-input-container label {
          color: #999999 !important;
        }

        #app_startup_loading_box {
          min-width: 320px;
          max-width: 500px;
          padding: 15px 20px;
          background: #d9edf7;
          border: 1px solid #bce8f1;
          border-radius: 4px;
          color: #31708f;
          font-size: 18px;
          box-shadow: 0 2px 8px rgba(0,0,0,0.15);
          text-align: center;
        }
        
        .filter-title-row {
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 10px;
          margin-top: 10px;
        }
        
        .filter-title-row h5 {
          margin: 0;
          font-size: 16px;
        }
        
        .filter-title-row .checkbox {
          margin: 0;
        }
        
        .filter-title-row .form-group {
          margin-bottom: 0;
        }

        #shiny-notification-panel {
          position: fixed;
          top: calc(50%);
          left: calc(50%);
          width: 500px;
          height: 200px;
          transform: translate(-50%, -50%);
          z-index: 9999;
        }

        .navbar-brand {
          display: flex !important;
          align-items: center !important;
          height: 60px !important;
        }

        .no-minor-ticks .irs-grid-pol.small {
          display: none;
        }

        .slider-right-fixed .irs-line {
          background: #428bca !important;
          border-color: #428bca !important;
        }

        .slider-right-fixed .irs-bar {
          background: #e5e5e5 !important;
          border-color: #e5e5e5 !important;
        }

        .slider-right-fixed .irs-bar-edge {
          display: none !important;
        }
        
        /* Sidebar section divider */
        .wide-sidebar-layout .sidebar-section-title {
          margin-top: 14px !important;
          margin-bottom: 12px !important;
          font-size: 20px !important;
          font-weight: 700 !important;
        }
        
        /* Make sidebar breaks compact but clear */
        .wide-sidebar-layout hr {
          margin-top: 18px !important;
          margin-bottom: 14px !important;
          border-top: 1px solid #cfcfcf !important;
        }
        
        /* Underline sidebar section title */
        .wide-sidebar-layout .sidebar-section-title {
          text-decoration: underline !important;
          text-underline-offset: 4px !important;
        }
        
        /* Align non-numeric filter checkboxes with the other filter titles */
        .wide-sidebar-layout #non-numeric_filters,
        .wide-sidebar-layout #patient_non-numeric_filters {
          padding-left: 8px !important;
          padding-right: 8px !important;
          box-sizing: border-box !important;
        }
        
        /* Make no-slider filter blocks align exactly like slider filter blocks */
        .wide-sidebar-layout #non-numeric_filters,
        .wide-sidebar-layout #patient_non-numeric_filters {
          padding-left: 0 !important;
          padding-right: 0 !important;
        }
        
        /* If old non-numeric checkbox styling still exists, neutralize it */
        .wide-sidebar-layout #non-numeric_filters > .checkbox,
        .wide-sidebar-layout #patient_non-numeric_filters > .checkbox {
          padding-left: 0 !important;
          margin-left: 0 !important;
        }
        
        /* Keep no-slider filter titles aligned with all other filter titles */
        .wide-sidebar-layout #non-numeric_filters .filter-block,
        .wide-sidebar-layout #patient_non-numeric_filters .filter-block {
          padding-left: 8px !important;
          padding-right: 8px !important;
          box-sizing: border-box !important;
        }
        
        /* Make their right-side checkbox text identical to other Enable checkboxes */
        .wide-sidebar-layout #non-numeric_filters .filter-title-row label,
        .wide-sidebar-layout #patient_non-numeric_filters .filter-title-row label {
          font-size: 15px !important;
          font-weight: normal !important;
          line-height: 1.3 !important;
          margin-bottom: 0 !important;
          white-space: nowrap !important;
        }
        
        /* Keep their titles identical to other filter titles */
        .wide-sidebar-layout #non-numeric_filters .filter-title-row h5,
        .wide-sidebar-layout #patient_non-numeric_filters .filter-title-row h5 {
          font-size: 17px !important;
          font-weight: 700 !important;
          line-height: 1.3 !important;
          margin: 0 !important;
        }
        
        /* Keep checkbox itself aligned nicely with bigger text */
        .wide-sidebar-layout #non-numeric_filters input[type='checkbox'],
        .wide-sidebar-layout #patient_non-numeric_filters input[type='checkbox'] {
          margin-top: 3px !important;
        }
        
        /* Compact filter blocks that only contain title + checkbox */
        .wide-sidebar-layout .filter-block .filter-control-area:empty {
          display: none !important;
        }
        
        .wide-sidebar-layout .filter-block:has(.filter-control-area:empty) {
          margin-top: 1px !important;
          margin-bottom: 30px !important;
        }
        
        /* Top sidebar section title */
        .wide-sidebar-layout .sidebar-section-title:first-child {
          margin-top: 0 !important;
        }
        

        .mode-tab {
          margin-top: 10px;
        }
      ")
    )
  ),
  
  tags$nav(
    class = "navbar navbar-default",
    style = "margin-bottom: 5px;",
    tags$div(
      class = "container-fluid",
      tags$div(
        class = "navbar-header",
        tags$span(
          class = "navbar-brand",
          style = "font-weight: bold; font-size: 30px;",
          HTML("ER<u>ASO</u>R")
        )
      )
    )
  ),
  
  div(
    id = "mode_tabs_wrapper",
    tabsetPanel(
      id = "mode_tabs",
      type = "tabs",
    
    ## =========================================================
    ## GENERAL KNOCKDOWN MODE
    ## =========================================================
    tabPanel(
      "General knockdown",
      div(
        class = "mode-tab",
        div(
          class = "wide-sidebar-layout",
          sidebarLayout(
            sidebarPanel(
            width = 3,
            fluidRow(
              column(
                12,
                
                h4(
                  "Input",
                  class = "sidebar-section-title"
                ),
                
                selectizeInput(
                  "ensemble_id_input",
                  label = tagList(
                    "Enter Gene or Ensembl ID (ENSG…):  ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "Please enter a gene, the script does not support transcript inputs.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  choices  = NULL,
                  multiple = FALSE,
                  options  = list(
                    placeholder = "Type a gene symbol (e.g. TP53) or Ensembl ID…",
                    maxOptions  = 50,
                    create      = TRUE
                  )
                ),
                
                hidden(
                checkboxInput(
                  "polymorphism_input",
                  label = tagList(
                    "Polymorphism analysis ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "This feature identifies single nucleotide polymorphisms (SNPs) for the selected Ensembl gene. SNPs are genetic variations in which one nucleotide has changed, for example, a C to a T.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  value = TRUE
                )),
                
                div(
                  # style = "display:none;", # this keeps the filter and button hidden
                  checkboxInput(
                    "linux_input",
                    label = tagList(
                      "Secondary structure calculation (Linux-OS only)",
                      tags$span(
                        tags$img(
                          src = "questionmark.png",
                          height = "20px",
                          style = "margin-bottom: 3px;"
                        ),
                        title = "Enable this setting only when running on a Linux operating system. Some features, such as the ViennaRNA analysis, require Linux. Enabling this option on other systems, such as Windows or Mac, may cause the program to fail.",
                        `data-placement` = "top",
                        `data-toggle` = "tooltip",
                        style = "cursor: pointer;"
                      )
                    ),
                    value = TRUE
                  )
                ),
                
                hr(),
                
                h4(
                  "Filter options",
                  class = "sidebar-section-title"
                ),
                
                div(
                  id = "non-numeric_filters",
                  
                  filterBlockUI(
                    title = "Filter out ASOs ending in G",
                    checkbox_id = "ASO_ending_G",
                    checkbox_value = TRUE,
                    tooltip = "Removes ASO sequences ending with G. These sequences have a higher chance to induce toxicity.",
                    content = NULL
                  ),
                  
                  filterBlockUI(
                    title = "Conserved in Mus musculus",
                    checkbox_id = "Conserved_input",
                    checkbox_value = FALSE,
                    tooltip = "This filter removes any ASOs where the target region is NOT conserved in the mouse genome. Enabling this filter will give ASOs that work in both human and mouse.",
                    content = NULL
                  )
                ),
                
                div(
                  id = "oligo_length_block",
                  class = "no-minor-ticks",
                  sliderInput(
                    "oligo_length_range",
                    label = tagList(
                      "ASO length: ",
                      tags$span(
                        tags$img(
                          src = "questionmark.png",
                          height = "20px",
                          style = "margin-bottom: 3px;"
                        ),
                        title = "Longer oligo lengths and a wider range of different lengths leads to longer runtime, we recommend a range of 3 lengths",
                        `data-placement` = "right",
                        `data-toggle` = "tooltip",
                        style = "cursor: pointer;"
                      )
                    ),
                    min = 15,
                    max = 25,
                    value = c(18, 19)
                  )
                ),
                
                div(
                  id = "filters",
                  fluidRow(
                    column(
                      12,
                      
                      filterBlockUI(
                        title = "GC content (%)",
                        checkbox_id = "gc_input",
                        checkbox_value = TRUE,
                        content = rangeFilterUI(
                          id    = "gc_content",
                          label = NULL,
                          min   = 0,
                          max   = 100,
                          value = c(40, 60),
                          step  = 1,
                          fixed = "none"
                        )
                      ),
                      
                      filterBlockUI(
                        title = "Acute neurotoxicity score (Hagedorn)",
                        checkbox_id = "tox_input",
                        checkbox_value = TRUE,
                        tooltip = "This quality control score estimates the potential toxicity of the complementary ASO targeting the mRNA sequence, providing insight into the safety risk of a chosen ASO. Higher values correspond to lower toxicity, which is favourable.",
                        content = rangeFilterUI(
                          id    = "tox_score",
                          label = NULL,
                          min   = 0,
                          max   = 136,
                          value = 60,
                          step  = 1,
                          fixed = "right",
                          label_right = "Minimum Tox. Score:"
                        )
                      ),
                      
                      filterBlockUI(
                        title = "Off-targets with perfect matches",
                        checkbox_id = "perfect_input",
                        checkbox_value = TRUE,
                        tooltip = "Number of off-targets, which are perfect complements to the ASO.",
                        content = div(
                          class = "no-minor-ticks",
                          rangeFilterUI(
                            id    = "perfect_hits",
                            label = NULL,
                            min   = 0,
                            max   = 5,
                            value = 1,
                            step  = 1,
                            fixed = "left",
                            label_left = "Maximum number of perfect matches:"
                          )
                        )
                      ),
                      
                
                      filterBlockUI(
                        title = "Off-targets with 1 mismatch",
                        checkbox_id = "mismatch_input",
                        checkbox_value = TRUE,
                        tooltip = "Desired number of off-targets with 1 mismatch/indel. Greater ranges lead to longer runtimes.",
                        content = rangeFilterUI(
                          id    = "mismatch_hits",
                          label = NULL,
                          min   = 0,
                          max   = 50,
                          value = 5,
                          step  = 1,
                          fixed = "left",
                          label_left = "Maximum number of off-targets with 1 mismatch:"
                        )
                      ),
                      
                      filterBlockUI(
                        title = "Polymorphism frequency",
                        checkbox_id = "Poly_input",
                        checkbox_value = TRUE,
                        tooltip = "This filter will remove any ASO sequences that fall outside of the determined window of polymorphism frequency. On the slider, only values that fall in the blue range are kept.",
                        content = rangeFilterUI(
                          id    = "pm_freq",
                          label = NULL,
                          min   = 0,
                          max   = 1,
                          value = c(0, 0.05),
                          step  = 0.01,
                          fixed = "none"
                        )
                      ),
                      
                      div(
                        style = "display:none;", # this keeps the filter hidden
                        h5(tagList(
                          HTML("<b>Accessibility </b>"),
                          tags$span(
                            tags$img(
                              src = "questionmark.png",
                              height = "20px",
                              style = "margin-bottom: 3px;"
                            ),
                            title = "This option is a quality control score for the accessibility of the target mRNA sequences. The score estimates how accessible the target mRNA sequences are, which is an important factor in finding potentially effective ASOs.",
                            `data-toggle` = "tooltip",
                            `data-placement` = "right",
                            style = "cursor: pointer;"
                          )
                        )),
                        # remove accessiblity as a filter btut keep it in the table
                        # fluidRow(
                        #   column(
                        #     9,
                        #     rangeFilterUI(
                        #       id    = "accessibility",
                        #       label = NULL,
                        #       min   = 0,
                        #       max   = 1,
                        #       value = c(0, 0.000001),
                        #       step  = 0.000001,
                        #       fixed = "none"
                        #     )
                        #   ),
                        #   column(3, checkboxInput("Accessibility_input", "Enable", value = FALSE))
                        # )
                      ),

                    )
                  )
                ),
                fluidRow(
                  column(4, actionButton("run_button", "Run")),
                  column(4, actionButton("reset_defaults", "Set to default"))
                )
              )
            )
          ),
          
          
          
          mainPanel(
            width = 9,
            tabsetPanel(
              id = "tabs_main_general",
              
              tabPanel(
                "Sequence results",
                div(
                  "This page displays the main results table, listing all ASOs that passed the filtering steps for the selected gene and providing information to help users select suitable ASO candidates.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                fluidRow(
                  column(
                    6,
                    p("This table shows the number of ASOs that did not meet filtering criterea."),
                    DTOutput("unfiltered_results_table")
                  ),
                  column(
                    2,
                    p("Download unfiltered results without off-target search here: "),
                    downloadButton("Download_unfiltered", "Download Unfiltered Results")
                  )
                ),
                hr(),
                fluidRow(
                  column(
                    3,
                    actionButton("toggle_cols", "Extended data")
                  ),
                  column(
                    4,
                    selectInput(
                      "main_table_sort_col",
                      "Sort by:",
                      choices = NULL
                    )
                  ),
                  column(
                    3,
                    selectInput(
                      "main_table_sort_dir",
                      "Direction:",
                      choices = c("Ascending" = "asc", "Descending" = "desc"),
                      selected = "asc"
                    )
                  )
                ),
                br(),
                p(tags$b(tags$u("Click")), " on a sequence to select it for RNaseH and off-target results. "),
                DT::dataTableOutput("results1"),
                hr(),
                downloadButton("Download_filtered", "Download Filtered Results")
              ),
              
              tabPanel(
                "RNase H cleavage results",
                div(
                  "This page predicts the optimal binding site for RNase H on the chosen target mRNA sequence. It uses a matrix of dinucleotide values and a sliding window approach to evaluate all possible binding sites along the mRNA sequence. The script accounts for chemical modifications on the ASO’s 5’ and 3’ ends, which can affect possible binding sites. Each potential site is scored based on the average from the dinucleotide matrix, and the results are ranked to highlight the most promising locations for RNase H cleavage.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                h3(textOutput("rnaseh_title")),
                div(uiOutput("rnaseh_info"), style = "margin-bottom: 15px;"),
                hr(),
                
                fluidRow(
                  column(
                    width = 4,
                    div(
                      style = "max-width: 320px;",
                      
                      numericInput(
                        "mod_5prime",
                        label = tagList(
                          "Amount of modified nucleotides at the 5' end  ",
                          tags$span(
                            tags$img(src = "questionmark.png", height = "20px"),
                            title = "This setting adds modifications to the 5' end of the ASO sequence, which affect the accessible binding sites on the target mRNA. An overlap of up to 2 modified nucleotides on the 3' end of the target mRNA is allowed.",
                            `data-toggle` = "tooltip",
                            style = "cursor: pointer;"
                          )
                        ),
                        value = 5,
                        min = 0,
                        max = 10,
                        width = "260px"
                      ),
                      
                      numericInput(
                        "mod_3prime",
                        label = tagList(
                          "Amount of modified nucleotides at the 3' end  ",
                          tags$span(
                            tags$img(src = "questionmark.png", height = "20px"),
                            title = "This setting adds modifications to the 3' end of the ASO sequence, which affect the accessible binding sites on the target mRNA. An overlap of up to 4 modified nucleotides on the 5' end of the target mRNA is allowed.",
                            `data-toggle` = "tooltip",
                            style = "cursor: pointer;"
                          )
                        ),
                        value = 5,
                        min = 0,
                        max = 10,
                        width = "260px"
                      ),
                      
                      div(
                        style = "margin-top: 10px;",
                        actionButton("add_mods", "Apply end modifications", class = "btn-primary")
                      ),
                      div(
                        style = "margin-top: 10px;",
                        downloadButton("download_rnaseh", "Download results")
                      )
                    )
                  ),
                  
                  column(
                    width = 8,
                    div(
                      style = "margin-top: -8px; display: flex; justify-content: flex-end;",
                      div(
                        style = "width: 620px; max-width: 100%; min-height: 170px; padding: 12px 18px; border: 1px solid #ddd; border-radius: 6px; background-color: #fafafa;",
                        h3("Visualised cleavage site:", style = "margin-top: 0; margin-bottom: 10px;"),
                        uiOutput("cleavage_visual")
                      )
                    )
                  )
                ),
                
                hr(),
                dataTableOutput("rnaseh_results")
              ),
              
              tabPanel(
                "Off target results",
                div(
                  "This page displays the off-targets for the selected target RNA sequence based on the allowed number of mismatches/indels. Off-target data is received via GGGenome. Protein Atlas is then to obtain  tissue expression data. This page also calculates the accessibility potential of each off-target.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                uiOutput("gggenome_status"),
                h3(textOutput("offtarget_title")),
                textOutput("aso_seq"),
                textOutput("numb_offtargets"),
                hr(),
                
                fluidRow(
                  column(
                    6,
                    selectInput(
                      "user_mismatch",
                      label = tagList(
                        "Select number of mismatches allowed  ",
                        tags$span(
                          tags$img(src = "questionmark.png", height = "20px"),
                          title = "This setting identifies off-target sequences according to the number of mismatches permitted in the target mRNA, with a default value of one mismatch. Allowing more mismatches results in more off-targets being detected but also increases the processing time.",
                          `data-toggle` = "tooltip",
                          style = "cursor: pointer;"
                        )
                      ),
                      choices = list(
                        "0" = 0,
                        "1" = 1,
                        "2" = 2,
                        "3" = 3
                      )
                    ),
                    actionButton("apply_mismatch", "Apply", class = "btn-primary")
                  ),
                  column(
                    6,
                    fluidRow("Run off-target tissue expression and OMIM disease search (may take some time)"),
                    fluidRow(
                      selectInput(
                        "target_tissue",
                        label = tagList(
                          "Select target tissue  ",
                          tags$span(
                            tags$img(src = "questionmark.png", height = "20px"),
                            title = "For each off-target, tissue-specific expression data is retrieved from the Protein Atlas. This data, combined with the off-target information, is then used to query OMIM for associated diseases. Currently, the analysis supports the brain, eyes, and liver, with additional tissues planned as further developments become available.",
                            `data-toggle` = "tooltip",
                            style = "cursor: pointer;"
                          )
                        ),
                        choices = c(
                          "Brain",
                          "Eye",
                          "Endrocrine tissue",
                          "Respiratory system",
                          "Proximal digestive tract",
                          "Gastrointestinal tract",
                          "Liver & galbladder",
                          "Pancreas",
                          "Kidney & urinary bladder",
                          "Male tissues",
                          "Female tissues",
                          "Muscle tissues",
                          "Connective & soft tissues",
                          "Skin",
                          "Bone marrow & lymphoid tissues"
                        )
                      ),
                      actionButton("PAtlas_OMIM_search", "Run")
                    )
                  )
                ),
                hr(),
                downloadButton("download_offtarget", "Download results", style = "margin-bottom: 15px;"),
                DTOutput("offtarget_results")
              )
            )
          )
        )
      )
    )
    ),
    
    ## =========================================================
    ## PATIENT-SPECIFIC MODE
    ## =========================================================
    
     if(show_patient_tab) tabPanel(
      "Patient-specific design",
      div(
        class = "mode-tab",
        div(
          class = "wide-sidebar-layout",
          sidebarLayout(
            sidebarPanel(
            width = 3,
            fluidRow(
              column(
                12,
                
                h4(
                  "Input",
                  class = "sidebar-section-title"
                ),
                
                selectizeInput(
                  "ensemble_id_input_patient",
                  label = tagList(
                    "Enter Gene or Ensembl ID (ENSG…):  ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "Please enter a gene, the script does not support transcript inputs.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  choices  = NULL,
                  multiple = FALSE,
                  options  = list(
                    placeholder = "Type a gene symbol (e.g. TP53) or Ensembl ID…",
                    maxOptions  = 50,
                    create      = TRUE
                  )
                ),
                
                selectInput(
                  "patient_analysis_mode",
                  label = tagList(
                    "Patient analysis mode ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "Choose whether to only inspect patient input and variants, generate both reference and SNV ASOs, or generate only patient-derived SNV ASOs.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  choices = c(
                    "Variant summary" = "summary_only",
                    "Variant summary + Variant ASOs" = "snv_only",
                    "Variant summary + Variant ASOs + General ASOs" = "summary_plus_aso"
                  ),
                  selected = "summary_only"
                ),
                
                fileInput(
                  "patient_variant_file",
                  label = tagList(
                    "Submit patient .bcf, .vcf, or .vcf.gz file: ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "ERASOR runs locally. Indexed .bcf and .vcf.gz files are best for later region-based reading. Plain .vcf is also accepted for now.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  accept = c(".bcf", ".vcf", ".vcf.gz")
                ),
                
                fileInput(
                  "patient_variant_index_file",
                  label = tagList(
                    "Upload optional index (.tbi or .csi): ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "Optional at this stage. Supported index formats are .tbi and .csi.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  accept = c(".tbi", ".csi")
                ),
                
                textInput(
                  "patient_region_input",
                  label = tagList(
                    "Optional chromosome region",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "Optional region in the format chr:start-end. If left empty, ERASOR will use the selected gene coordinates plus the flank below.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  placeholder = "e.g. chr11:532242-537343"
                ),
                
                numericInput(
                  "patient_flank",
                  "Flank around selected gene (bp)",
                  value = 0,
                  min = 0,
                  max = 1000000000,
                  step = 1
                ),
                
                div(
                  id = "patient_optional_analysis_controls",
                  
                  checkboxInput(
                    "patient_filter_snv_asos",
                    label = "Apply filters to patient-derived SNV ASOs",
                    value = TRUE
                  ),
                  
                  hidden(
                  checkboxInput(
                    "patient_polymorphism_input",
                    label = tagList(
                      "Polymorphism analysis ",
                      tags$span(
                        tags$img(
                          src = "questionmark.png",
                          height = "20px",
                          style = "margin-bottom: 3px;"
                        ),
                        title = "This feature identifies single nucleotide polymorphisms (SNPs) for the selected Ensembl gene. SNPs are genetic variations in which one nucleotide has changed, for example, a C to a T.",
                        `data-placement` = "right",
                        `data-toggle` = "tooltip",
                        style = "cursor: pointer;"
                      )
                    ),
                    value = TRUE
                  )),
                  
                  checkboxInput(
                    "patient_linux_input",
                    label = tagList(
                      "Secondary structure calculation (Linux-OS only)",
                      tags$span(
                        tags$img(
                          src = "questionmark.png",
                          height = "20px",
                          style = "margin-bottom: 3px;"
                        ),
                        title = "Enable this setting only when running on a Linux operating system. Some features, such as the ViennaRNA analysis, require Linux. Enabling this option on other systems may cause the program to fail.",
                        `data-placement` = "top",
                        `data-toggle` = "tooltip",
                        style = "cursor: pointer;"
                      )
                    ),
                    value = TRUE
                  )
                ),
                
                div(
                  id = "patient_filter_options_container",
                
                hr(),
                
                h4(
                  "Filter options",
                  class = "sidebar-section-title"
                ),
                
                div(
                  id = "patient_non-numeric_filters",
                  
                  
                  filterBlockUI(
                    title = "Filter out ASOs ending in G",
                    checkbox_id = "patient_ASO_ending_G",
                    checkbox_value = TRUE,
                    tooltip = "Removes ASO sequences ending with G. These sequences have a higher chance to induce toxicity.",
                    content = NULL
                  ),
                  
                  filterBlockUI(
                    title = "Conserved in Mus musculus",
                    checkbox_id = "patient_Conserved_input",
                    checkbox_value = FALSE,
                    tooltip = "This filter removes any ASOs where the target region is NOT conserved in the mouse genome. Enabling this filter will give ASOs that work in both human and mouse. ",
                    content = NULL
                  )
                  
                ),
                
                div(
                  id = "patient_oligo_length_block",
                  class = "no-minor-ticks",
                  sliderInput(
                    "patient_oligo_length_range",
                    label = tagList(
                      "ASO length: ",
                      tags$span(
                        tags$img(
                          src = "questionmark.png",
                          height = "20px",
                          style = "margin-bottom: 3px;"
                        ),
                        title = "Longer oligo lengths and a wider range of different lengths leads to longer runtime, we recommend a range of 3 lengths",
                        `data-placement` = "right",
                        `data-toggle` = "tooltip",
                        style = "cursor: pointer;"
                      )
                    ),
                    min = 15,
                    max = 25,
                    value = c(18, 19)
                  )
                ),
                
                filterBlockUI(
                  title = "GC content (%)",
                  checkbox_id = "patient_gc_input",
                  checkbox_value = TRUE,
                  content = rangeFilterUI(
                    id    = "patient_gc_content",
                    label = NULL,
                    min   = 0,
                    max   = 100,
                    value = c(40, 60),
                    step  = 1,
                    fixed = "none"
                  )
                ),
                
                filterBlockUI(
                  title = "Acute neurotoxicity score (Hagedorn)",
                  checkbox_id = "patient_tox_input",
                  checkbox_value = TRUE,
                  tooltip = "This quality control score estimates the potential toxicity of the complementary ASO targeting the mRNA sequence. Higher values correspond to lower toxicity, which is favourable.",
                  content = rangeFilterUI(
                    id    = "patient_tox_score",
                    label = NULL,
                    min   = 0,
                    max   = 136,
                    value = 60,
                    step  = 1,
                    fixed = "right",
                    label_right = "Minimum Tox. Score:"
                  )
                ),
                
                filterBlockUI(
                  title = "Off-targets with perfect matches",
                  checkbox_id = "patient_perfect_input",
                  checkbox_value = TRUE,
                  tooltip = "Number of off-targets, which are perfect complements to the ASO.",
                  content = div(
                    class = "no-minor-ticks",
                    rangeFilterUI(
                      id    = "patient_perfect_hits",
                      label = NULL,
                      min   = 0,
                      max   = 5,
                      value = 1,
                      step  = 1,
                      fixed = "left",
                      label_left = "Maximum number of perfect matches:"
                    )
                  )
                ),
                
                filterBlockUI(
                  title = "Off-targets with 1 mismatch",
                  checkbox_id = "patient_mismatch_input",
                  checkbox_value = TRUE,
                  tooltip = "Desired number of off-targets with 1 mismatch/indel. Greater ranges lead to longer runtimes.",
                  content = rangeFilterUI(
                    id    = "patient_mismatch_hits",
                    label = NULL,
                    min   = 0,
                    max   = 50,
                    value = 5,
                    step  = 1,
                    fixed = "left",
                    label_left = "Maximum number of off-targets with 1 mismatch:"
                  )
                ),
                
                filterBlockUI(
                  title = "Polymorphism frequency",
                  checkbox_id = "patient_Poly_input",
                  checkbox_value = TRUE,
                  tooltip = "This filter will remove any ASO sequences that fall outside of the determined window of polymorphism frequency. On the slider, only values that fall in the blue range are kept.",
                  content = rangeFilterUI(
                    id    = "patient_pm_freq",
                    label = NULL,
                    min   = 0,
                    max   = 1,
                    value = c(0, 0.05),
                    step  = 0.01,
                    fixed = "none"
                  )
                )
                ),
                
                fluidRow(
                  column(4, actionButton("run_button_patient", "Run")),
                  column(4, actionButton("reset_defaults_patient", "Set to default"))
                )
              )
            )
          ),
          
          mainPanel(
            width = 9,
            tabsetPanel(
              id = "tabs_main_patient",
              
              tabPanel(
                "Variant summary",
                div(
                  "This page checks the patient variants and their positions in the selected gene.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                h4("Patient variant summary"),
                tableOutput("patient_summary"),
                hr(),
                downloadButton("Download_unfiltered_patient", "Download Variant Summary"),
                hr(),
                h4("Detected SNVs in selected region"),
                DTOutput("patient_snvs_table")
              ),
              
              tabPanel(
                "Sequence results",
                div(
                  "This page shows reference-gene ASOs and, separately, patient-derived SNV ASOs.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                h4("Reference-gene ASO results"),
                p("These ASOs are generated from the selected reference gene sequence."),
                
                fluidRow(
                  column(
                    3,
                    actionButton("toggle_cols_patient", "Extended data")
                  ),
                  column(
                    4,
                    selectInput(
                      "patient_main_table_sort_col",
                      "Sort by:",
                      choices = NULL
                    )
                  ),
                  column(
                    3,
                    selectInput(
                      "patient_main_table_sort_dir",
                      "Direction:",
                      choices = c("Ascending" = "asc", "Descending" = "desc"),
                      selected = "asc"
                    )
                  )
                ),
                br(),
                
                DTOutput("patient_results1"),
                br(),
                downloadButton(
                  "Download_unfiltered_patient_reference",
                  "Download Unfiltered Reference ASO Results"
                ),
                downloadButton(
                  "Download_filtered_patient_reference",
                  "Download Reference ASO Results"
                ),
                hr(),
                h4("Patient-derived SNV ASO results"),
                p("These ASOs are generated from windows overlapping a single patient SNV. Metadata columns are shown only in this table."),
                
                fluidRow(
                  column(
                    3,
                    actionButton("toggle_cols_patient_snv", "Extended data")
                  ),
                  column(
                    4,
                    selectInput(
                      "patient_snv_table_sort_col",
                      "Sort by:",
                      choices = NULL
                    )
                  ),
                  column(
                    3,
                    selectInput(
                      "patient_snv_table_sort_dir",
                      "Direction:",
                      choices = c("Ascending" = "asc", "Descending" = "desc"),
                      selected = "asc"
                    )
                  )
                ),
                br(),
                
                DTOutput("patient_snv_results"),
                br(),
                downloadButton("Download_patient_snv_results", "Download Patient SNV ASO Results"),
                
                hr(),
                h4("Ambiguous unphased multi-variant ASO windows"),
                DTOutput("patient_ambiguous_results"),
                br(),
                downloadButton(
                  "Download_patient_ambiguous_results",
                  "Download Ambiguous Windows"
                )
              ),
              
              tabPanel(
                "RNase H cleavage results",
                div(
                  "This page predicts the optimal RNase H binding site for the selected patient-tab ASO.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                h3(textOutput("patient_rnaseh_title")),
                div(uiOutput("patient_rnaseh_info"), style = "margin-bottom: 15px;"),
                hr(),
                
                fluidRow(
                  column(
                    width = 4,
                    div(
                      style = "max-width: 320px;",
                      
                      numericInput(
                        "patient_mod_5prime",
                        "Amount of modified nucleotides at the 5' end",
                        value = 5,
                        min = 0,
                        max = 10,
                        width = "260px"
                      ),
                      
                      numericInput(
                        "patient_mod_3prime",
                        "Amount of modified nucleotides at the 3' end",
                        value = 5,
                        min = 0,
                        max = 10,
                        width = "260px"
                      ),
                      
                      div(
                        style = "margin-top: 10px;",
                        actionButton(
                          "patient_add_mods",
                          "Apply end modifications",
                          class = "btn-primary"
                        )
                      ),
                      
                      div(
                        style = "margin-top: 10px;",
                        downloadButton(
                          "patient_download_rnaseh",
                          "Download results"
                        )
                      )
                    )
                  ),
                  
                  column(
                    width = 8,
                    div(
                      style = "margin-top: -8px; display: flex; justify-content: flex-end;",
                      div(
                        style = "width: 620px; max-width: 100%; min-height: 170px; padding: 12px 18px; border: 1px solid #ddd; border-radius: 6px; background-color: #fafafa;",
                        h3("Visualised cleavage site:", style = "margin-top: 0; margin-bottom: 10px;"),
                        uiOutput("patient_cleavage_visual")
                      )
                    )
                  )
                ),
                
                hr(),
                DTOutput("patient_rnaseh_results")
              ),
              
              tabPanel(
                "Off target results",
                div(
                  "This page displays the off-targets for the selected target RNA sequence based on the allowed number of mismatches/indels. Off-target data is received via GGGenome. Protein Atlas is then to obtain  tissue expression data. This page also calculates the accessibility potential of each off-target.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                h3(textOutput("patient_offtarget_title")),
                textOutput("patient_aso_seq"),
                textOutput("patient_numb_offtargets"),
                hr(),
                
                fluidRow(
                  column(
                    6,
                    selectInput(
                      "patient_user_mismatch",
                      "Select number of mismatches allowed",
                      choices = list(
                        "0" = 0,
                        "1" = 1,
                        "2" = 2,
                        "3" = 3
                      ),
                      selected = 2
                    ),
                    actionButton(
                      "patient_apply_mismatch",
                      "Apply",
                      class = "btn-primary"
                    )
                  ),
                  
                  column(
                    6,
                    fluidRow("Run off-target tissue expression and OMIM disease search."),
                    fluidRow(
                      selectInput(
                        "patient_target_tissue",
                        "Select target tissue",
                        choices = c(
                          "Brain",
                          "Eye",
                          "Endrocrine tissue",
                          "Respiratory system",
                          "Proximal digestive tract",
                          "Gastrointestinal tract",
                          "Liver & galbladder",
                          "Pancreas",
                          "Kidney & urinary bladder",
                          "Male tissues",
                          "Female tissues",
                          "Muscle tissues",
                          "Connective & soft tissues",
                          "Skin",
                          "Bone marrow & lymphoid tissues"
                        )
                      ),
                      actionButton(
                        "patient_PAtlas_OMIM_search",
                        "Run"
                      )
                    )
                  )
                ),
                
                hr(),
                downloadButton(
                  "patient_download_offtarget",
                  "Download results",
                  style = "margin-bottom: 15px;"
                ),
                DTOutput("patient_offtarget_results")
              )
              
            )
          )
        )
      )
    )
    ),
    
    
    ## =========================================================
    ## SPECIFIC ASO INPUT MODE
    ## =========================================================
    tabPanel(
      "Specific ASO input",
      div(
        class = "mode-tab",
        div(
          class = "wide-sidebar-layout",
          sidebarLayout(
            sidebarPanel(
            width = 3,
            
            div(
              style = "margin: 15px 0; padding: 15px; border: 2px solid #5bc0de; border-radius: 8px; background-color: #f4fbfd;",
              
              h4("Specific ASO input"),
              
              p("Analyse manually entered ASO sequences without generating all ASOs from a reference gene, with a reference gene, or against patient-specific consensus sequences."),
              
              selectInput(
                "single_aso_analysis_mode",
                label = tagList(
                  "Specific ASO analysis mode ",
                  tags$span(
                    tags$img(
                      src = "questionmark.png",
                      height = "20px",
                      style = "margin-bottom: 3px;"
                    ),
                    title = "Choose whether submitted ASOs should be analysed without a reference gene, with full reference-gene annotations, or against patient-specific consensus sequences for conservation and polymorphism frequency.",
                    `data-placement` = "right",
                    `data-toggle` = "tooltip",
                    style = "cursor: pointer;"
                  )
                ),
                choices = c(
                  "ASOs without known target" = "no_reference",
                  "ASOs targeting a reference gene sequence" = "with_reference_full",
                  "ASOs targeting a gene with patient variants" = "with_reference_consensus_cf"
                ),
                selected = "no_reference"
              ),
              checkboxInput(
                "single_aso_linux_input",
                label = tagList(
                  "Secondary structure calculation (Linux-OS only)",
                  tags$span(
                    tags$img(
                      src = "questionmark.png",
                      height = "20px",
                      style = "margin-bottom: 3px;"
                    ),
                    title = "Enable this setting only when running on a Linux operating system. This calculates ASO self-folding energy and ASO duplex energy using ViennaRNA. Enabling this option on systems without ViennaRNA may cause the analysis to fail.",
                    `data-placement` = "top",
                    `data-toggle` = "tooltip",
                    style = "cursor: pointer;"
                  )
                ),
                value = TRUE
              ),
              
              conditionalPanel(
                condition = "input.single_aso_analysis_mode != 'no_reference'",
                
                selectizeInput(
                  "ensemble_id_input_single_aso",
                  label = tagList(
                    "Reference Gene or Ensembl ID (ENSG…): ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "For full annotation, submitted ASOs are mapped to this selected reference gene. For patient consensus mode, this gene defines the reference region and coordinate system.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  choices = NULL,
                  multiple = FALSE,
                  options = list(
                    placeholder = "Type a gene symbol or Ensembl ID…",
                    maxOptions = 50,
                    create = TRUE
                  )
                ),
                
                hidden(
                checkboxInput(
                  "single_aso_polymorphism_input",
                  "Polymorphism analysis",
                  value = TRUE
                ))
              ),
              
              conditionalPanel(
                condition = "input.single_aso_analysis_mode == 'with_reference_consensus_cf'",
                
                fileInput(
                  "single_aso_variant_file",
                  label = tagList(
                    "Submit patient .bcf, .vcf, or .vcf.gz file: ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "Used to build patient-specific reference, haplotype, and ambiguous consensus sequences for ASO matching.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  accept = c(".bcf", ".vcf", ".vcf.gz")
                ),
                
                fileInput(
                  "single_aso_variant_index_file",
                  "Upload optional index (.tbi or .csi):",
                  accept = c(".tbi", ".csi")
                ),
                
                numericInput(
                  "single_aso_consensus_max_ambiguous",
                  "Maximum ambiguous consensus possibilities",
                  value = 32,
                  min = 1,
                  max = 512,
                  step = 1
                )
              ),
              
              textAreaInput(
                "aso_seq_input",
                label = "Enter ASO sequences",
                value = "",
                placeholder = "Enter ASO sequences separated by comma, point, newline, or whitespace.",
                rows = 8
              ),
              
              tags$small("Only sequences containing A, C, G, and T are recognized."),
              br(), br(),
              
              strong("Detected ASO sequences:"),
              verbatimTextOutput("parsed_asos"),
              
              br(),
              
              actionButton(
                "run_button_single_aso",
                "Run specific ASO analysis",
                class = "btn-primary"
              )
            )
          ),
          
          mainPanel(
            width = 9,
            tabsetPanel(
              id = "tabs_main_single_aso",
              
              tabPanel(
                "Sequence results",
                div(
                  "This page displays the output table for manually submitted ASO sequences.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                
                fluidRow(
                  column(
                    3,
                    actionButton("toggle_cols_single_aso", "Extended data")
                  ),
                  column(
                    4,
                    selectInput(
                      "single_aso_table_sort_col",
                      "Sort by:",
                      choices = NULL
                    )
                  ),
                  column(
                    3,
                    selectInput(
                      "single_aso_table_sort_dir",
                      "Direction:",
                      choices = c("Ascending" = "asc", "Descending" = "desc"),
                      selected = "asc"
                    )
                  )
                ),
                
                br(),
                
                p("Click on a sequence to select it for RNase H and off-target results."),
                DT::dataTableOutput("single_aso_results"),
                hr(),
                downloadButton("Download_single_aso_results", "Download Results")
              ),
              
              tabPanel(
                "RNase H cleavage results",
                div(
                  "This page predicts the optimal RNase H binding site for the selected manually submitted ASO.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                
                h3(textOutput("single_aso_rnaseh_title")),
                div(uiOutput("single_aso_rnaseh_info"), style = "margin-bottom: 15px;"),
                hr(),
                
                fluidRow(
                  column(
                    width = 4,
                    div(
                      style = "max-width: 320px;",
                      
                      numericInput(
                        "single_aso_mod_5prime",
                        "Amount of modified nucleotides at the 5' end",
                        value = 5,
                        min = 0,
                        max = 10,
                        width = "260px"
                      ),
                      
                      numericInput(
                        "single_aso_mod_3prime",
                        "Amount of modified nucleotides at the 3' end",
                        value = 5,
                        min = 0,
                        max = 10,
                        width = "260px"
                      ),
                      
                      div(
                        style = "margin-top: 10px;",
                        actionButton(
                          "single_aso_add_mods",
                          "Apply end modifications",
                          class = "btn-primary"
                        )
                      ),
                      
                      div(
                        style = "margin-top: 10px;",
                        downloadButton(
                          "single_aso_download_rnaseh",
                          "Download results"
                        )
                      )
                    )
                  ),
                  
                  column(
                    width = 8,
                    div(
                      style = "margin-top: -8px; display: flex; justify-content: flex-end;",
                      div(
                        style = "width: 620px; max-width: 100%; min-height: 170px; padding: 12px 18px; border: 1px solid #ddd; border-radius: 6px; background-color: #fafafa;",
                        h3("Visualised cleavage site:", style = "margin-top: 0; margin-bottom: 10px;"),
                        uiOutput("single_aso_cleavage_visual")
                      )
                    )
                  )
                ),
                
                hr(),
                DTOutput("single_aso_rnaseh_results")
              ),
              
              tabPanel(
                "Off target results",
                div(
                  "This page displays the off-targets for the selected target RNA sequence based on the allowed number of mismatches/indels. Off-target data is received via GGGenome. Protein Atlas is then to obtain  tissue expression data. This page also calculates the accessibility potential of each off-target.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                
                h3(textOutput("single_aso_offtarget_title")),
                textOutput("single_aso_seq"),
                textOutput("single_aso_numb_offtargets"),
                hr(),
                
                fluidRow(
                  column(
                    6,
                    selectInput(
                      "single_aso_user_mismatch",
                      "Select number of mismatches allowed",
                      choices = list(
                        "0" = 0,
                        "1" = 1,
                        "2" = 2,
                        "3" = 3
                      ),
                      selected = 2
                    ),
                    actionButton(
                      "single_aso_apply_mismatch",
                      "Apply",
                      class = "btn-primary"
                    )
                  ),
                  
                  column(
                    6,
                    fluidRow("Run off-target tissue expression and OMIM disease search."),
                    fluidRow(
                      selectInput(
                        "single_aso_target_tissue",
                        "Select target tissue",
                        choices = c(
                          "Brain",
                          "Eye",
                          "Endrocrine tissue",
                          "Respiratory system",
                          "Proximal digestive tract",
                          "Gastrointestinal tract",
                          "Liver & galbladder",
                          "Pancreas",
                          "Kidney & urinary bladder",
                          "Male tissues",
                          "Female tissues",
                          "Muscle tissues",
                          "Connective & soft tissues",
                          "Skin",
                          "Bone marrow & lymphoid tissues"
                        )
                      ),
                      actionButton(
                        "single_aso_PAtlas_OMIM_search",
                        "Run"
                      )
                    )
                  )
                ),
                
                hr(),
                downloadButton(
                  "single_aso_download_offtarget",
                  "Download results",
                  style = "margin-bottom: 15px;"
                ),
                DTOutput("single_aso_offtarget_results")
              )
            )
          )
          )
        )
      )
    ),
    
    ## =========================================================
    ## SPECIFIC ASO INPUT MODE
    ## =========================================================
    tabPanel(
      "Help & Citation",
      div(
        class = "mode-tab",)
    )
    
    
    )
  ),
  
  tags$script(HTML("
  var loadingDotsInterval = setInterval(function() {
    var el = document.getElementById('loading_dots');
    if (!el) return;

    if (el.textContent === '.') {
      el.textContent = '..';
    } else if (el.textContent === '..') {
      el.textContent = '...';
    } else {
      el.textContent = '.';
    }
  }, 500);

  setTimeout(function() {

    var txt = document.getElementById('startup_loading_text');
    var dots = document.getElementById('loading_dots');
    var box = document.getElementById('app_startup_loading_box');

    /* stop the dot animation */
    clearInterval(loadingDotsInterval);

    if (txt) {
      txt.textContent = 'GGGenome is unavailable. Please try again another time.';
    }

    if (dots) {
      dots.textContent = '';
    }

    if (box) {
      box.style.background = '#f2dede';
      box.style.border = '1px solid #ebccd1';
      box.style.color = '#a94442';
    }

  }, 120000);
  
  function markCompletedUploads() {
    $('.shiny-input-container input[type=file]').each(function() {
      var container = $(this).closest('.shiny-input-container');
      var bar = container.find('.progress-bar');

      if (bar.length) {
        var width = bar.attr('style') || '';
        if (width.indexOf('100%') !== -1) {
          container.addClass('file-upload-complete');
        } else {
          container.removeClass('file-upload-complete');
        }
      }
    });
  }

  $(document).on('change', 'input[type=file]', function() {
    var container = $(this).closest('.shiny-input-container');
    container.removeClass('file-upload-complete');
  });

  setInterval(markCompletedUploads, 300);
")),
  
  tags$script(
    HTML('$(function () { $("[data-toggle=\'tooltip\']").tooltip(); });')
  )
)
