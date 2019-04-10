def top_variables(loadings, pc, n=20):
    return loadings[pc].abs().nlargest(n)
