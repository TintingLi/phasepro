


# place the code and data to run only once -------------------------------
library(shiny)


library(magrittr)
library(data.table)
library(Biostrings)
library(RColorBrewer)
palette(c(brewer.pal(8, 'Set2')))

source('function.R')

#fa sequence
#fa <- readAAStringSet('~/Desktop/projects/phase_separation/results/fa/human/all_protein_coding_genes_without_asterix_aa_sequence.fa')
#saveRDS(fa, file = './data/all_pro_aas_without_asterix.rda')
fa <- readRDS('data/all_pro_aas_without_asterix.rda')

#low complex domain scores and AA sequences:
lcd_list <- readRDS("data/lcd_list.rda")
#motif locus:
motifs <- readRDS('data/motifs_clean.rda')
#annotations:
pro_anno <- readRDS('data/proteins_anno.rda')



# Define UI for app that draws a histogram ----
ui <- fluidPage(
	
	# App title ----
	titlePanel(
		list(HTML('<a href="http://www.nki.nl/"><img src="ncba.jpg" class="img.unframed" align="right" alt="NKI" height="auto" width="15%"></a>'), "phasepro: phase separation predictioin"),
			   windowTitle="phasepro: phase separation prediction"),
	
	# Sidebar layout with input and output definitions ----
	sidebarLayout(
		
		# Sidebar panel for inputs ----
		sidebarPanel(
			
			textInput(inputId = 'gene_symbol',
					  label = 'HGNC gene symbol:',
					  value = 'TP53'),
			
		submitButton(text = 'update view')

		),
		
		# Main panel for displaying outputs ----
		mainPanel(
			
			tabsetPanel(
				
				tabPanel('introduction',
					uiOutput('introduction')),
				tabPanel('plot',
						 plotOutput(outputId = "Plot", height = '900px') ),
				tabPanel('sequence and downloads',
					h3('data download:'),
					downloadButton('downloadfa', 'protein sequence'),
					br(),
					br(),
					downloadButton('downloadps', 'LCD scores'),
					
					h3('protein sequence:'),
					verbatimTextOutput('protein_seqs'),
					h3('LCD scores:'),
					tableOutput('protein_lcd'),
					br()),
				
				tabPanel('public links',
				  uiOutput('links')
				  )
	
			)
			

		)
	)
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
	
	#hgnc gene symbols:
	i = reactive({ input$gene_symbol })
	
	#seqs and LCD scores:
	output$protein_seqs <- renderPrint(	fa[names(fa) == i()] %>% as.character)
	output$protein_lcd <- renderTable(lcd_list[[i()]])
	output$downloadfa <- downloadHandler(
		filename = function(){
			paste0(i(), '.fa')
		},
		content = function(file){
			writeXStringSet(fa[names(fa) == i()], file)
		}
	)

	output$downloadps <- downloadHandler(
		filename = function(){
			paste0(i(), '_lcd_score.csv')
		},
		content = function(file){
			fwrite(lcd_list[[i()]], file)
		}
	)
	
	
	#public links:
	output$links <- renderUI({
		tagList(
			h3(paste('HGNC gene symbol:', i())),
			p(a('GeneCards', 
			  href = pro_anno[gene_symbol == i(), gene_cards_website])),

			p(a('Uniprot', 
			  href = pro_anno[gene_symbol == i(), uniprot_website])),

			p(a('Google', 
			  href = pro_anno[gene_symbol == i(), google_website])),

			p(a('Google Scholar', 
			  href = pro_anno[gene_symbol == i(), google_scholar_website])),

			p(a('Google image', 
			  href = pro_anno[gene_symbol == i(), google_image_website]))
		)
	})

	#introduction html:
	output$introduction <- renderUI({
		tagList(
			h3('Purpose'),
			p('This is a web tool facilities predicting protein potential of phase separation.'),
			br(),
			
			h3('Instructions'),
			
			p('Enter the HGNC gene symbol and click the "update view" button.'),
			h4('plot window'),
			
			p(),
			h4('sequence and downloads window'),
			p(),
			h4('public links window'),
			p(),
			
			h3('Softwares'),
			p()
		)
	})
	
		
	#main plot:
	output$Plot <- renderPlot({
		
		plot_LCD_Motif_AAs(motif.dt = motifs[gene_symbol == i()],
						   lcd.dt = lcd_list[[i()]], 
						   gene_symbol = i())
		
	})
	
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)


