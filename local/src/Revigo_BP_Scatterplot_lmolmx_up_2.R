# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000280","nuclear division",1.81594680672738,4.13018662865131,-5.27710532532376,2.51321760006794,-15.357456539712,0.840927165114428,0.0035574),
c("GO:0001523","retinoid metabolic process",0.480527462703246,-2.26875014933033,6.53253546491148,1.93951925261862,-1.6233809067082,0.880132074909494,0.48245368),
c("GO:0006260","DNA replication",1.10633066994468,2.30612999440054,5.60796841348347,2.29885307640971,-4.18568094333325,0.92591036974005,0.20391537),
c("GO:0006325","chromatin organization",3.4083924680114,4.08598960971421,-3.63840689828437,2.78604121024255,-1.65563140864112,0.937410689205305,0.22001979),
c("GO:0006833","water transport",0.139688215902106,-6.16099928418685,-2.33061659204264,1.41497334797082,-1.65039089804791,0.966322922795289,0.18582266),
c("GO:0007059","chromosome segregation",1.63714589037269,6.52932388503017,-0.227145926781756,2.46834733041216,-15.357456539712,0.72857461967453,0),
c("GO:0007586","digestion",0.6034530926971,-0.754255748141711,-7.82117742562276,2.03742649794062,-1.78439391492121,0.995525535440476,0),
c("GO:0008203","cholesterol metabolic process",0.659328379057943,-0.490030212389814,6.23320252193929,2.07554696139253,-3.40180044478591,0.801668470345558,0),
c("GO:0008544","epidermis development",1.76007152036654,5.7635814410151,2.67228108090208,2.4996870826184,-2.88669069280976,0.980911462105051,0),
c("GO:0009812","flavonoid metabolic process",0.0726378722690954,-5.85253157563644,-5.6987022200451,1.14612803567824,-3.0100610841237,0.975937068474964,0.05751192),
c("GO:0009913","epidermal cell differentiation",1.1007431413086,5.39812359141013,2.51668797744782,2.29666519026153,-2.50922756971621,0.970033573401961,0.49530939),
c("GO:0015701","bicarbonate transport",0.167625859082528,-6.89113195994461,-1.61409898557042,1.49136169383427,-1.85729759822865,0.95188492600778,0.18582266),
c("GO:0015849","organic acid transport",1.55333296083142,-6.33723940809922,-0.827003700602121,2.4456042032736,-1.37539677757811,0.952368409688803,0.37358773),
c("GO:0016266","O-glycan processing",0.240263731351623,3.19424357404257,5.53489270594379,1.64345267648619,-1.70308519957066,0.91285098711341,0.32035131),
c("GO:0016441","post-transcriptional gene silencing",0.217913616807286,-2.71624786842165,3.64411091515665,1.60205999132796,-3.72025549607623,0.918797672908859,0.10604009),
c("GO:0019216","regulation of lipid metabolic process",1.86623456445214,0.144221258307355,1.38772366958391,2.52504480703685,-1.85729759822865,0.957295201591366,0.16137261),
c("GO:0019218","regulation of steroid metabolic process",0.514052634519752,1.92006679737993,1.00827932844064,1.96848294855394,-4.85374459177175,0.926990742464779,0.03036141),
c("GO:0021700","developmental maturation",1.50304520310667,5.40072820040314,3.71977419679181,2.43136376415899,-1.37539677757811,0.991490378414642,0.19007907),
c("GO:0031055","chromatin remodeling at centromere",0.0447002290886741,5.8864101373323,-3.14950769873069,0.954242509439325,-1.45586100916277,0.957402553490369,0.13866768),
c("GO:0034032","purine nucleoside bisphosphate metabolic process",0.68167849360228,0.948220914496649,6.04503608783424,2.0899051114394,-1.70308519957066,0.82618605562996,0.43197127),
c("GO:0035383","thioester metabolic process",0.508465105883668,1.06009482383905,-7.72741927809098,1.96378782734556,-1.85729759822865,0.919823795713771,0.06837386),
c("GO:0042044","fluid transport",0.189975973626865,-6.92296690700225,0.357594018513858,1.54406804435028,-1.56719825593107,0.965509707869087,0.19033689),
c("GO:0042634","regulation of hair cycle",0.156450801810359,-4.35246606738765,-6.54215090896018,1.46239799789896,-2.48506921401123,0.981961246542,0.02688298),
c("GO:0042908","xenobiotic transport",0.268201374532044,-7.279591478763,-0.587165576403982,1.69019608002851,-1.3346317494155,0.964552745333635,0.19767729),
c("GO:0045787","positive regulation of cell cycle",1.97798513717383,-1.58966763312341,-4.93136247505991,2.55022835305509,-8.88672607156992,0.724766652641527,0),
c("GO:0046717","acid secretion",0.139688215902106,-5.10937563469907,-1.15631323804819,1.41497334797082,-1.93646575563909,0.966322922795289,0),
c("GO:0048285","organelle fission",1.95563502262949,4.41309126160325,-4.85373347413212,2.54530711646582,-15.357456539712,0.899067364432866,0.37233798),
c("GO:0050000","chromosome localization",0.458177348158909,-5.39662790440332,0.584729837690648,1.91907809237607,-1.87743471760161,0.945488633284445,0.16984106),
c("GO:0051276","chromosome organization",2.43616248533274,3.7389852744495,-4.76658816310368,2.64048143697042,-4.68779843484001,0.896818317498234,0.38664003),
c("GO:0051302","regulation of cell division",1.05604291221993,-2.76631756802284,-7.56045920609681,2.27875360095283,-2.46095888481802,0.974327411517944,0.0413688),
c("GO:0051304","chromosome separation",0.0614628149969269,6.81879250233727,-0.340470331591302,1.07918124604762,-7.02157224147822,0.771400577308842,0.49356138),
c("GO:0051783","regulation of nuclear division",0.810191652232218,-5.3000630874143,5.28939738462259,2.16435285578444,-6.7585959404082,0.884106543002424,0.04011424),
c("GO:0051988","regulation of attachment of spindle microtubules to kinetochore",0.100575515449517,-1.30661121279949,-4.60838157989429,1.27875360095283,-2.48506921401123,0.756473408976317,0.47942313),
c("GO:0052695","cellular glucuronidation",0.111750572721685,0.0293558248287265,6.54016563033193,1.32221929473392,-1.65039089804791,0.908901809807186,0.36587675),
c("GO:0071103","DNA conformation change",0.508465105883668,3.44692115053124,-5.52148195977628,1.96378782734556,-5.70702639051271,0.870218117152348,0.31834856),
c("GO:0071459","protein localization to chromosome, centromeric region",0.139688215902106,-6.65607152076524,1.71529719685918,1.41497334797082,-1.85729759822865,0.954956306196119,0.16984106),
c("GO:1902930","regulation of alcohol biosynthetic process",0.268201374532044,0.171329023672369,-1.09724339021779,1.69019608002851,-3.11988504058249,0.949146572877424,0.10800743),
c("GO:2001251","negative regulation of chromosome organization",0.441414762250657,-4.8223905058515,5.17053265141474,1.90308998699194,-1.68973594059807,0.815814183463218,0.4908475));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


chosen <- read_excel("/scratch/trcanmed/biobanca/local/share/data/chosen_revigo_signature.xlsx", 
                     sheet = "up")
colnames(chosen) <- "description"
chosen$chosen <- "YES"
#chosen <- rbind(chosen , c("muscle organ development", "YES"))

one.data <- merge(one.data, chosen, by = "description", all.x = TRUE)
one.data$chosen[is.na(one.data$chosen)] <- "NO"
one.data <- one.data[, c(2,1,3,4,5,6,7,8,9,10)]
rownames(one.data) <- one.data$term_ID
one.data$new_log_size <- NA
one.data["GO:0015849", "chosen"] <- "YES"
one.data["GO:0007059", "chosen"] <- "YES"
one.data["GO:0000280", "chosen"] <- "YES"
	
### impongo anzichè la mediana il primo quartile e al posto di 2.5 il terzo quartile

for (i in rownames(one.data)) {
  if (one.data[i, "log_size"] < 1.96) {
    one.data[i, "new_log_size"] <- 3
  } else if (one.data[i, "log_size"] > 1.96 & one.data[i, "log_size"] < 2.3 | one.data[i, "log_size"] == 1.96) {
    one.data[i, "new_log_size"] <- 5
  } else {
    one.data[i, "new_log_size"] <- 8
  }
}

# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below
textSize <- 5
largerSize <- textSize + 2


unmute_theme <- theme_bw() +
  theme(
    text = element_text(size = textSize, family='sans'),
    axis.title = element_text(size = largerSize),
    axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
    axis.text.y = element_text(size = textSize, color="black"),
    plot.title = element_text(size = largerSize, hjust = 0.5),
    legend.title = element_text(size=largerSize),
    legend.text = element_text(size=textSize)
  )


p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point(aes( plot_X, plot_Y, colour = value, size = new_log_size))# + scale_size_area();
p1 <- p1 + scale_colour_gradient(low = "red", high = "yellow")#, limits = c( min(one.data$value), 5) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = new_log_size), shape = 21, fill = "transparent", colour = I (alpha ("white", 0.6) ))# + scale_size_area();
#p1 <- p1 + scale_size(one.data$new_log_size) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
p1 <- p1 + scale_size_continuous(labels = c(1.96, 2.3, 3.0), breaks = c(3, 5, 8)) + theme_bw()
ex <- one.data [ one.data$chosen == "YES", ];
#ex <- one.data [ one.data$dispensability < 0.15, ]
p1 <- p1 + geom_label_repel(data = ex, aes(plot_X, plot_Y, label = description), color = 'black', size = (6/2.8),
                            box.padding = unit(0.35, "lines"),
                            point.padding = unit(0.5, "lines"),
                            segment.color = 'grey50'
)+unmute_theme;
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y", size = "log_size");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
p1 <- p1 + unmute_theme
ggsave("revigo_lmo_lmx_up_def.pdf", width = 120, height = 107, useDingbats=FALSE, units = "mm")
p1
dev.off()
