library(RColorBrewer)
library(ggplot2)
nice_theme = list(theme_bw(), theme(legend.position='bottom'))

patient_colors <- list(
    'HIV status'=brewer.pal(2, "Set3")[1:2],
    'Tuberculosis status'=rev(brewer.pal(4, "Greens")),
    'Meningitis'=ggthemes::tableau_color_pal('Tableau 10')(5),
    'Is RNA outlier'=c('grey95', 'red'),
    'Cluster'=brewer.pal(5, "Accent")[1:5],
    'Sex'=c('lightblue', 'darkgoldenrod2'),
    'Steroids treatment'=c('green', 'grey90'),
    'TB treatment'=c('green', 'grey90', 'grey96'),
    'CSF RNA batch'=c(brewer.pal(2, "Set2")[1:2], 'grey96')
)

patient_colors_values <- list(
    'HIV status'=c('Positive', 'Negative'),
    'Tuberculosis status'=c('Definite', 'Probable', 'Possible', '-'),
    # TODO: rename Tuberculosis to Tuberculous
    'Meningitis'=c('Tuberculosis', 'Cryptococcal', 'Viral', 'Healthy control', 'Bacterial'),
    'Is RNA outlier'=c('False', 'True'),
    'Cluster'=c('1', '2', '3', '4', '5'),
    'Sex'=c('M', 'F'),
    'Steroids treatment'=c('Yes', 'No'),
    'TB treatment'=c('Yes', 'No', '?'),
    'CSF RNA batch'=c('FC1', 'FC2',  '-')
)

for(group in names(patient_colors)) {
    names(patient_colors[[group]]) <- patient_colors_values[[group]]
}



color_tuberculosis = list(
    scale_color_manual(values=patient_colors$`Tuberculosis status`, name='Tuberculosis status')
)

color_meningitis = list(
    scale_color_manual(values=patient_colors$Meningitis, name='Meningitis')
)
fill_meningitis = list(
    scale_fill_manual(values=patient_colors$Meningitis, name='Meningitis')
)