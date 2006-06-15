from rasmus import fff
reload(fff)


findex = fff.FeatureIndex()
findex.read("motifs2.fff")

f = findex.getFeatures("yeast", "chromosome4", 1e3, 2e3)
#l = findex.lookup("yeast", "chromosome4", 1e3, 2e3)


