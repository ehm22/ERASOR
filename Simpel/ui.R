library(DT)
library(shiny)
library(shinyjs)
library(shinythemes)

ui <- fluidPage(
  useShinyjs(),
  theme = shinythemes::shinytheme("flatly"),
  tags$head(tags$style(
    HTML(
      "
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
    "
    )
  )),
  tags$nav(
    class = "navbar navbar-default",
    style = "margin-bottom: 5px;",
    tags$div(class = "container-fluid", tags$div(
      class = "navbar-header",
      tags$span(class = "navbar-brand", style = "font-weight: bold; font-size: 30px;", HTML("ER<u>ASO</u>R"))
    ))
  ),
  sidebarLayout(
    sidebarPanel(fluidRow(
      column(
        12,
        textInput(
          "ensemble_id_input",
          label = tagList(
            "Enter Ensembl id:  ",
            tags$span(
              tags$img(
                src = "questionmark.png",
                height = "20px",
                style = "margin-bottom: 3px;"
              ),
              title = "Please enter a valid Ensembl Gene ID (ENSG…) for the gene of interest.",
              `data-placement` = "right",
              `data-toggle` = "tooltip",
              style = "cursor: pointer;"
            )
          ),
          value = "ENSG00000100284"
        ),
        checkboxInput(
          "single_aso_input",
          label = "Single ASO analysis, enter ASO sequence:",
          value = FALSE
        ),
        conditionalPanel(
          condition = "input.single_aso_input == true",
          textInput(
            "aso_seq_input",
            label = "",
            value = "",
            placeholder = "e.g., GCCTCAGTCTGCTTCGCACC"
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
              title = "Filters ASO sequences ending with G. These sequences have a higher chance to induce toxicity. ",
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
            "Conserved in mouse ",
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
          value = TRUE
        ),
        sliderInput(
          "oligo_length_range",
          label = tagList(
            "Oligo length: ",
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
        ),
        div(
          id = "filters",
        fluidRow(
          column(
            12,
            
            h5(tagList(
              HTML("<b>Off-targets with perfect matches</b> "),
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = "Desired number of perfectly complementary off-targets. ",
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            )),
            fluidRow(
              column(
                4,
                selectInput(
                  "dropdown_input_a",
                  "",
                  selected = "==",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_a", "", value = 1, min = 0)),
              checkboxInput("perfect_input", "Enable", value = TRUE)
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
                4,
                selectInput(
                  "dropdown_input_b",
                  "",
                  selected = "<=",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_b", "", value = 10, min = 0)),
              checkboxInput("mismatch_input", "Enable", value = TRUE)
            ),
            
            h5(tagList(
              HTML("<b>Accessibility </b>"),
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = "This option is a quality control score for the accessibility of the target mRNA sequences. The score estimates how accessible the target mRNA sequences are, which is an important factor in finding potentially effective ASOs. ",
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            )),
            fluidRow(
              column(
                4,
                selectInput(
                  "dropdown_input_c",
                  "",
                  selected = ">",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_c", "", value = 1E-6, min = 0)),
              checkboxInput("Accessibility_input", "Enable", value =
                              TRUE)
            ),
            
            h5(tagList(
              HTML("<b>Polymorphism frequence </b>"),
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = "This quality control score estimates the probability that a polymorphism (SNP) is present in the target mRNA sequences. Lower values are preferred, as they indicate a reduced likelihood of polymorphic variation.",
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            )),
            fluidRow(
              column(
                4,
                selectInput(
                  "dropdown_input_d",
                  "",
                  selected = "<",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_d", "", value = 0.05, min = 0, step = 0.01)),
              checkboxInput("Poly_input", "Enable", value = TRUE)
            ),
            
            h5(tagList(
              HTML("<b>Toxicity score </b>"),
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = "This quality control score estimates the potential toxicity of the complementary ASO targeting the mRNA sequence, providing insight into the safety risk of a chosen ASO. Higher values correspond to lower toxicity, which is favourable.",
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            )),
            fluidRow(
              column(
                4,
                selectInput(
                  "dropdown_input_e",
                  "",
                  selected = ">=",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_e", "", value = 60, min = 0)),
              checkboxInput("tox_input", "Enable", value = TRUE)
            ),
            fluidRow(
              column(6, actionButton("run_button", "Run")),
              column(6, actionButton("reset_defaults", "Set to default", style = "margin-left: -33px;"))
            )),
          )
        )
      )
    ),width = 3),
    
    # Updated main panel for RNase H script.
    mainPanel(tabsetPanel(
      id = "tabs_main",
      tabPanel(
        "Sequence results",
        div(
          "This is the main page of the application, which displays the primary output table. The table lists all potential target mRNA sequences for the provided Ensembl ID and includes information to help the user select suitable ASO targets. After selecting an ASO, the application automatically redirects to the RNase H tab, and the chosen ASO becomes available for further analysis on all other tabs.",
          style = "margin-top: 20px; font-size: 18px;"
        ),
        hr(),
        fluidRow(
          column(6,
                 p("This table shows the number of ASOs that did not meet filtering criterea."),
                 DTOutput("unfiltered_results_table")
                ),
          column(2,
                 p("Download unfiltered results without off-target search here: "),
                 downloadButton("Download_unfiltered", "Download Unfiltered Results")
               )),
        hr(),
        fluidRow(
          column(
            3,
            actionButton("toggle_cols", "Extended data")
          )
        ),
        br(),
        DT::dataTableOutput('results1'),
        hr(),
        downloadButton("Download_filtered", "Download Filtered Results"),
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
        
        fluidRow(column(
          6,
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
            value = 0,
            min = 0,
            max = 10,
            width = "60%"
          ),
        ), column(
          6,
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
            value = 0,
            min = 0,
            max = 10,
            width = "60%"
          )
        )),
        actionButton("add_mods", "Apply end modifications", class = "btn-primary"),
        hr(),
        
        downloadButton("download_rnaseh", "Download results", style = "margin-bottom: 15px;"),
        dataTableOutput("rnaseh_results"),
        hr(),
        
        h3("Visualised cleavage site: "),
        uiOutput("cleavage_visual"),
      ),
      
      tabPanel(
        "Off target results", 
        div("This page displays the off-targets for the selected target mRNA sequence based on the allowed number of mismatches. For the complementary ASO sequence, it shows mismatches, deletions, insertions, protein name and additional information from GGGenome. The protein name is then given to Protein Atlas to retrieve the tissue expression data. This page also calculates the accessibility potential of each off-target.", style="margin-top: 20px; font-size: 18px;"),
        hr(),
        uiOutput("gggenome_status"),
        h3(textOutput("offtarget_title")),
        textOutput("aso_seq"),
        textOutput("numb_offtargets"),
        hr(),
        fluidRow(
          column(6, 
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
                 actionButton("apply_mismatch", "Apply", class = "btn-primary"),
          ),
          column(6,
                 fluidRow("Run off-target tissue expression and OMIM disease search (may take some time)"),
                 
                 fluidRow(
                   selectInput("target_tissue", label = tagList(
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
                   actionButton("PAtlas_OMIM_search", "Run"))
          )
        ),
        hr(),
        downloadButton("download_offtarget", "Download results", style = "margin-bottom: 15px;"),
        DTOutput("offtarget_results"),
      )
    ),
    width = 9
  )
  ),
  tags$script(
    HTML(
      '$(function () { $("[data-toggle=\'tooltip\']").tooltip(); });'
    )
  )
)