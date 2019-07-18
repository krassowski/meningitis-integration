library(RColorBrewer)
library(ggplot2)
nice_theme = list(theme_bw(), theme(legend.position='bottom'))

patient_colors <- list(
    'HIV status'=brewer.pal(2, "Set3")[1:2],
    'Tuberculosis status'=rev(brewer.pal(4, "Greens")),
    'Meningitis'=ggthemes::tableau_color_pal('Tableau 10')(5),
    'Is RNA outlier'=c('grey95', 'red')

)

patient_colors_values <- list(
    'HIV status'=c('Positive', 'Negative'),
    'Tuberculosis status'=c('Definite', 'Probable', 'Possible', '-'),
    # TODO: rename Tuberculosis to Tuberculous
    'Meningitis'=c('Tuberculosis', 'Cryptococcal', 'Viral', 'Healthy control', 'Bacterial'),
    'Is RNA outlier'=c('False', 'True')
)

for(group in names(patient_colors)) {
    names(patient_colors[[group]]) <- patient_colors_values[[group]]
}



color_meningitis = list(
    scale_color_manual(values=patient_colors$Meningitis, name='Meningitis')
)
fill_meningitis = list(
    scale_fill_manual(values=patient_colors$Meningitis, name='Meningitis')
)