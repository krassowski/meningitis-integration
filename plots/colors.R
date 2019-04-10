library(RColorBrewer)

patient_colors <- list(
    'HIV status'=brewer.pal(2, "Set3")[1:2],
    'Tuberculosis status'=rev(brewer.pal(4, "Greens")),
    'Meningitis'=ggthemes::tableau_color_pal('Tableau 10')(4)
)

patient_colors_values <- list(
    'HIV status'=c('Positive', 'Negative'),
    'Tuberculosis status'=c('Definite', 'Probable', 'Possible', '-'),
    'Meningitis'=c('Tuberculosis', 'Cryptococcal', 'Viral', 'Healthy control')
)

for(group in names(patient_colors)) {
    names(patient_colors[[group]]) <- patient_colors_values[[group]]
}
