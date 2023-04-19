#CARREGANDO PACOTES NECESSÁRIOS (INSTALAR SE FOR PRECISO)
library(readxl)
library(tidyr)
library(dplyr)
library(tidyverse)
library(pheatmap)


## MANIPULAÇÃO DOS DADOS E PREPARAÇÃO PARA A ANÁLISE
########
# Carrega o sample dos dados e renomeaia as colunas do dataframe para "genome" e "homologous_group"
df <- read_excel("COLAR_AQUI_O_DIRETÓRIO_ATUAL/sample_data_acineto.xlsx")

df <- rename(df, genome = orgsymbol, homologous_group = ko)


# Criaa uma coluna com valor "1" para contar a presença do homologous_group em cada genoma analisado
df <- df %>%
  mutate(count = 1L) # alterando o tipo da coluna "count" para integer

# Caso ocorra algum valor repetido, esse comando exclui para não interferir na analise
df <- unique(df)

# Cria uma nova tabela utilizando a função pivot_wider() do pacote tidyr atribuindo valor 0 para genes ausentes entre as bactérias
df_wide <- df %>%
  pivot_wider(names_from = genome, values_from = count, values_fill = 0)


# Remove a coluna de homologous_group e renomeaa as colunas de genome para incluir o prefixo "genome:"
df_wide <- df_wide %>%
  select(-homologous_group) %>%
  rename_all(~paste0("genome:", .))

########


##PIPELINE PARA ANÁLISE DO GENE CONTENT
########

# Transforma os dados do genoma numa matriz
gene_groups <- as.matrix(df_wide)

# Normaliza a matriz dividindo cada valor pela soma de todos os valores 
normalized_matrix <- gene_groups / sum(gene_groups)

# Criando a matriz de correlação
cor_matrix <- cor(normalized_matrix, method="pearson")

colnames(cor_matrix) <- gsub("genome:", "", colnames(cor_matrix))
rownames(cor_matrix) <- gsub("genome:", "", rownames(cor_matrix))



# Plotando o heatmap
pheatmap(cor_matrix,
         #Parametros que alteram tamanho, angulo e coloca a fonte em negrito
         fontsize_row = 15, fontsize_col = 15,
         fontfamily_row = "bold", fontfamily_col = "bold",
         angle_col = 45)
########
