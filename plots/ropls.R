library(ropls)
plot = ropls::plot

plot_summary = function(opls_result) {
    par(mfrow=c(2,2))
    plot(opls_result, typeVc="overview", parDevNewL=F)
    plot(opls_result, typeVc="outlier", parDevNewL=F)
    plot(opls_result, typeVc="x-score", parDevNewL=F)
    plot(opls_result, typeVc="x-loading", parDevNewL=F)
}

plot_opls_all = function(opls_result, test=F) {
    par(mfrow=c(3,3))
    plot(opls_result, typeVc="overview", parDevNewL=F)
    plot(opls_result, typeVc="outlier", parDevNewL=F)
    plot(opls_result, typeVc="x-score", parDevNewL=F)
    plot(opls_result, typeVc="x-loading", parDevNewL=F)
    plot(opls_result, typeVc="xy-score", parDevNewL=F)
    plot(opls_result, typeVc="xy-weight", parDevNewL=F) # w is grey, C is black
    plot(opls_result, typeVc="x-variance", parDevNewL=F)
    plot(opls_result, typeVc="correlation", parDevNewL=F)
    #plot(opls_result, typeVc="permutation", parDevNewL=F, permI=2)

    # for DA
    # plot(opls_result, typeVc="predict-train", parDevNewL=F)
        
    if (test)
        plot(opls_result, typeVc="predict-test", parDevNewL=F)
}


plot_opls = function(opls_result, typeVc, ...) {
    plot(opls_result, typeVc=typeVc, parDevNewL=F, ...)
}

# workaround to the very original idea of manipulating device and window in plot() functions
run_and_plot_opls = function(..., fun=plot){
    options(device="png", mar=c(0,0,0,0), dpi=400, width=500, res=1500)
    result = fun(...)
    dev.off()
    img = png::readPNG(paste(getwd(), 'Rplot001.png', sep='/'))
    par(bg=NA,mar=c(0,0,0,0),oma=c(0,0,0,0))
    plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
    rasterImage(img,0,0,1,1, interpolate=F)
    result
}

# TODO: bi-plot