library(DT)
library(shiny)
library(shinyjs)
library(shinythemes)
##&*## just to troll annelot
library(bslib)
##&*## just to troll annelot

# if this is TRUE the patient specific tab will show, if FALSE then the tab is hidden
show_patient_tab <- FALSE 

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

ui <- fluidPage(
  useShinyjs(),
  
  div(
    id = "app_startup_loading",
    div(
      id = "app_startup_loading_box",
      span(id = "startup_loading_text", "Please wait. Loading"),
      span(id = "loading_dots", ".")
    )
  ),
  
    theme = shinythemes::shinytheme("flatly"),
  
  #### troll purple for annelot ##&*####
  # theme = bs_theme(
  #   bootswatch = "flatly",
  #   primary = "purple",
  #   bg = "#ffe4ec",   # light pink background
  #   fg = "#2b002b"    # dark text so it stays readable
  # ),
  #### troll purple for annelot ##&*####
  
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
        
        /* rows clickable */
          table.dataTable tbody tr {
          cursor: pointer;
          transition: background-color 0.15s ease;
        }

        /* hover highlight */
          table.dataTable tbody tr:hover {
          background-color: #e8f4ff !important;
          }
        
        /* center DT column headers */
            table.dataTable thead th {
           text-align: center !important;
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
        sidebarLayout(
          sidebarPanel(
            width = 3,
            fluidRow(
              column(
                12,
                
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
                
                checkboxInput(
                  "single_aso_input",
                  label = "Get characteristics for specific ASOs",
                  value = FALSE
                ),
                
                conditionalPanel(
                  condition = "input.single_aso_input == true",
                  tagList(
                    textAreaInput(
                      "aso_seq_input",
                      label = NULL,
                      value = "",
                      placeholder = "Enter ASO sequences separated by comma, point or white space.",
                      rows = 5
                    ),
                    tags$small(
                      "Only sequences containing A, C, G, and T are recognized. "
                    ),
                    br(),
                    strong("Detected ASO sequences:"),
                    verbatimTextOutput("parsed_asos")
                  )
                ),
                
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
                ),
                
                checkboxInput(
                  "ASO_ending_G",
                  label = tagList(
                    "Filter out ASOs ending in G ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "Removes ASO sequences ending with G. These sequences have a higher chance to induce toxicity.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  value = TRUE
                ),
                
                checkboxInput(
                  "Conserved_input",
                  label = tagList(
                    "Conserved in Mus musculus ",
                    tags$span(
                      tags$img(
                        src = "questionmark.png",
                        height = "20px",
                        style = "margin-bottom: 3px;"
                      ),
                      title = "This feature identifies conserved regions and orthologous genes in the mouse genome. Enabling this filter will show ASOs which could be used in mouse models.",
                      `data-placement` = "right",
                      `data-toggle` = "tooltip",
                      style = "cursor: pointer;"
                    )
                  ),
                  value = FALSE
                ),
                div(
                  style = "display:none;", # this keeps the filter and button hidden
                  checkboxInput(
                    "linux_input",
                    label = tagList(
                      "Accessibility calculation (Linux-OS only)",
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
                    value = FALSE
                  )
                ),
                
                div(
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
                    value = c(18, 20)
                  )
                ),
                
                div(
                  id = "filters",
                  fluidRow(
                    column(
                      12,
                      
                      h5(tagList(
                        HTML("<b>GC content (%) </b>")
                      )),
                      fluidRow(
                        column(
                          9,
                          rangeFilterUI(
                            id    = "gc_content",
                            label = NULL,
                            min   = 0,
                            max   = 100,
                            value = c(40, 60),
                            step  = 1,
                            fixed = "none"
                          )
                        ),
                        column(
                          3,
                          checkboxInput("gc_input", "Enable", value = TRUE)
                        )
                      ),
                      
                      h5(tagList(
                        HTML("<b>Acute neurotoxicity score (Hagedorn) </b>"),
                        tags$span(
                          tags$img(
                            src = "questionmark.png",
                            height = "20px",
                            style = "margin-bottom: 3px;"
                          ),
                          title = "This quality control score estimates the potential toxicity of the complementary ASO targeting the mRNA sequence, providing insight into the safety risk of a chosen ASO. Higher values correspond to lower toxicity, which is favourable.",
                          `data-toggle` = "tooltip",
                          style = "cursor: pointer;"
                        )
                      )),
                      fluidRow(
                        column(
                          9,
                          rangeFilterUI(
                            id    = "tox_score",
                            label = NULL,
                            min   = 0,
                            max   = 136,
                            value = 60,
                            step  = 1,
                            fixed = "right",
                            label_right = "Minimum Tox. Score"
                          )
                        ),
                        column(
                          3,
                          checkboxInput("tox_input", "Enable", value = TRUE)
                        )
                      ),
                      
                      h5(tagList(
                        HTML("<b>Off-targets with perfect matches</b> "),
                        tags$span(
                          tags$img(
                            src = "questionmark.png",
                            height = "20px",
                            style = "margin-bottom: 3px;"
                          ),
                          title = "Number of off-targets, which are perfect complements to the ASO.",
                          `data-toggle` = "tooltip",
                          `data-placement` = "right",
                          style = "cursor: pointer;"
                        )
                      )),
                      fluidRow(
                        div(
                          class = "no-minor-ticks",
                          column(
                            9,
                            rangeFilterUI(
                              id    = "perfect_hits",
                              label = NULL,
                              min   = 0,
                              max   = 5,
                              value = 1,
                              step  = 1,
                              fixed = "left",
                              label_left = "Maximum number of perfect matches"
                            )
                          ),
                          column(3, checkboxInput("perfect_input", "Enable", value = TRUE))
                        )
                      ),
                      
                      h5(tagList(
                        HTML("<b>Off-targets with 1 mismatch </b>"),
                        tags$span(
                          tags$img(
                            src = "questionmark.png",
                            height = "20px",
                            style = "margin-bottom: 3px;"
                          ),
                          title = "Desired number of off-targets with 1 mismatch/indel. Greater ranges lead to longer runtimes.",
                          `data-toggle` = "tooltip",
                          `data-placement` = "right",
                          style = "cursor: pointer;"
                        )
                      )),
                      fluidRow(
                        column(
                          9,
                          rangeFilterUI(
                            id    = "mismatch_hits",
                            label = NULL,
                            min   = 0,
                            max   = 50,
                            value = 5,
                            step  = 1,
                            fixed = "left",
                            label_left = "Maximum number of off-targets with 1 mismatch"
                          )
                        ),
                        column(3, checkboxInput("mismatch_input", "Enable", value = TRUE))
                      ),
                      
                      h5(tagList(
                        HTML("<b>Polymorphism frequency </b>"),
                        tags$span(
                          tags$img(
                            src = "questionmark.png",
                            height = "20px",
                            style = "margin-bottom: 3px;"
                          ),
                          title = "This quality control score estimates the probability that a polymorphism (SNP) is present in the target mRNA sequences. Lower values are preferred.",
                          `data-toggle` = "tooltip",
                          `data-placement` = "right",
                          style = "cursor: pointer;"
                        )
                      )),
                      fluidRow(
                        column(
                          9,
                          rangeFilterUI(
                            id    = "pm_freq",
                            label = NULL,
                            min   = 0,
                            max   = 1,
                            value = 0.05,
                            step  = 0.01,
                            fixed = "left"
                          )
                        ),
                        column(3, checkboxInput("Poly_input", "Enable", value = TRUE))
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
                        fluidRow(
                          column(
                            9,
                            rangeFilterUI(
                              id    = "accessibility",
                              label = NULL,
                              min   = 0,
                              max   = 1,
                              value = c(0, 0.000001),
                              step  = 0.000001,
                              fixed = "none"
                            )
                          ),
                          column(3, checkboxInput("Accessibility_input", "Enable", value = FALSE))
                        )
                      ),
                      
                      fluidRow(
                        column(4, actionButton("run_button", "Run")),
                        column(4, actionButton("reset_defaults", "Set to default"))
                      )
                    )
                  )
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
                  "This is the main page of the application, which displays the primary output table. The table lists all potential target mRNA sequences for the provided Ensembl ID and includes information to help the user select suitable ASO targets. After selecting an ASO, the application automatically redirects to the RNase H tab, and the chosen ASO becomes available for further analysis on all other tabs.",
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
                  )
                ),
                br(),
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
                  "This page displays the off-targets for the selected target mRNA sequence based on the allowed number of mismatches. For the complementary ASO sequence, it shows mismatches, deletions, insertions, protein name and additional information from GGGenome. The protein name is then given to Protein Atlas to retrieve the tissue expression data. This page also calculates the accessibility potential of each off-target.",
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
    ),
    
    ## =========================================================
    ## PATIENT-SPECIFIC MODE
    ## =========================================================
    
     if(show_patient_tab) tabPanel(
      "Patient-specific design",
      div(
        class = "mode-tab",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            fluidRow(
              column(
                12,
                
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
                  max = 5000,
                  step = 1
                ),
                
                fluidRow(
                  column(
                    6,
                    actionButton("run_button_patient", "Run")
                  )
                )
              )
            )
          ),
          
          mainPanel(
            width = 9,
            tabsetPanel(
              id = "tabs_main_patient",
              
              tabPanel(
                "Input summary",
                div(
                  "This page currently only checks the patient input setup. It verifies the selected gene, stages the uploaded variant and index files, resolves the chromosome region, and shows a summary of the detected input. No variant reading or ASO generation is performed yet.",
                  style = "margin-top: 20px; font-size: 18px;"
                ),
                hr(),
                h4("Patient input summary"),
                tableOutput("patient_summary"),
                hr(),
                h4("Status checks"),
                DTOutput("unfiltered_results_table_patient"),
                hr(),
                h4("Resolved inputs"),
                DTOutput("results1_patient"),
                hr(),
                downloadButton("Download_unfiltered_patient", "Download Input Summary"),
                br(), br(),
                downloadButton("Download_filtered_patient", "Download Input Summary"),
                hr(),
                h4("Detected SNVs in selected region"),
                DTOutput("patient_snvs_table"),
                
                hr(),
                h4("Reference sequence"),
                verbatimTextOutput("patient_reference_sequence"),
                
                hr(),
                h4("Patient-specific sequence"),
                verbatimTextOutput("patient_modified_sequence")
              )
            )
          )
        )
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

  }, 60000);
")),
  
  tags$script(
    HTML('$(function () { $("[data-toggle=\'tooltip\']").tooltip(); });')
  )
)
