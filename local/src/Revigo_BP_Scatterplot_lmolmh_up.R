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
revigo.data <- rbind(c("GO:0000280","nuclear division",1.81594680672738,2.87891853359783,6.04091524977036,2.51321760006794,-15.357456539712,0.840927165114428,0.0035574),
c("GO:0006260","DNA replication",1.10633066994468,-1.12801233388299,-6.5728822281773,2.29885307640971,-4.18568094333325,0.92591036974005,0.18517545),
c("GO:0006325","chromatin organization",3.4083924680114,3.22195385596304,4.48883253788291,2.78604121024255,-1.65563140864112,0.937410689205305,0.22001979),
c("GO:0006833","water transport",0.139688215902106,-5.0014822672987,-0.37489360907018,1.41497334797082,-1.65039089804791,0.966322922795289,0.18582266),
c("GO:0007059","chromosome segregation",1.63714589037269,5.93149317776697,-2.87443957203104,2.46834733041216,-15.357456539712,0.72857461967453,0),
c("GO:0007586","digestion",0.6034530926971,-4.23875478658737,6.75490642297986,2.03742649794062,-1.78439391492121,0.995525535440476,0),
c("GO:0008202","steroid metabolic process",1.43040733083757,-1.57238321227466,-5.06106936816062,2.40993312333129,-4.52986151029481,0.881817738835022,0.48843231),
c("GO:0008544","epidermis development",1.76007152036654,3.50537078872693,-4.9187416837839,2.4996870826184,-2.88669069280976,0.980911462105051,0),
c("GO:0009812","flavonoid metabolic process",0.0726378722690954,-6.83525659417024,4.50011434907719,1.14612803567824,-3.0100610841237,0.975937068474964,0.0533159),
c("GO:0009913","epidermal cell differentiation",1.1007431413086,3.85989430635333,-5.07821752774142,2.29666519026153,-2.50922756971621,0.970033573401961,0.49530939),
c("GO:0015701","bicarbonate transport",0.167625859082528,-6.00291100741668,0.47378730985469,1.49136169383427,-1.85729759822865,0.95188492600778,0.18582266),
c("GO:0015849","organic acid transport",1.55333296083142,-6.28485323399009,-0.458867269344027,2.4456042032736,-1.37539677757811,0.952368409688803,0.37358773),
c("GO:0016126","sterol biosynthetic process",0.251438788623792,-0.881388045998378,-5.3809209133945,1.66275783168157,-6.41181374698838,0.829690630874989,0),
c("GO:0016266","O-glycan processing",0.240263731351623,-0.17268031247715,-6.90093231433511,1.64345267648619,-1.70308519957066,0.91285098711341,0.32035131),
c("GO:0016441","post-transcriptional gene silencing",0.217913616807286,3.17618328917196,-0.699982084823821,1.60205999132796,-3.72025549607623,0.918797672908859,0.10604009),
c("GO:0019216","regulation of lipid metabolic process",1.86623456445214,6.61360333944735,0.211579485618287,2.52504480703685,-1.85729759822865,0.957295201591366,0.16137261),
c("GO:0019218","regulation of steroid metabolic process",0.514052634519752,5.55431551896015,0.532777002884661,1.96848294855394,-4.85374459177175,0.926990742464779,0.03036141),
c("GO:0021700","developmental maturation",1.50304520310667,3.05290794274066,-5.94022889365493,2.43136376415899,-1.37539677757811,0.991490378414642,0.19007907),
c("GO:0031055","chromatin remodeling at centromere",0.0447002290886741,5.05093947228045,4.34749402059179,0.954242509439325,-1.45586100916277,0.957402553490369,0.13866768),
c("GO:0034032","purine nucleoside bisphosphate metabolic process",0.68167849360228,-2.80615047097355,-6.15548163603748,2.0899051114394,-1.70308519957066,0.82618605562996,0.22718674),
c("GO:0035383","thioester metabolic process",0.508465105883668,-0.73834123836661,7.75327228355359,1.96378782734556,-1.85729759822865,0.919823795713771,0.06252381),
c("GO:0042044","fluid transport",0.189975973626865,-7.08110932032198,-1.34549226835952,1.54406804435028,-1.56719825593107,0.965509707869087,0.19033689),
c("GO:0042634","regulation of hair cycle",0.156450801810359,-2.47398770577952,7.35796440040483,1.46239799789896,-2.48506921401123,0.981961246542,0.02688298),
c("GO:0042908","xenobiotic transport",0.268201374532044,-7.30937152636051,-0.214434626771579,1.69019608002851,-1.3346317494155,0.964552745333635,0.19767729),
c("GO:0045787","positive regulation of cell cycle",1.97798513717383,-2.01771166262418,3.83124195320692,2.55022835305509,-8.88672607156992,0.724766652641527,0),
c("GO:0046717","acid secretion",0.139688215902106,-6.98714004465448,1.018103114525,1.41497334797082,-1.93646575563909,0.966322922795289,0),
c("GO:0048285","organelle fission",1.95563502262949,2.48510055435717,5.47176895387326,2.54530711646582,-15.357456539712,0.899067364432866,0.37233798),
c("GO:0050000","chromosome localization",0.458177348158909,-5.59638259014728,-1.95405788432327,1.91907809237607,-1.87743471760161,0.945488633284445,0.16984106),
c("GO:0051276","chromosome organization",2.43616248533274,3.16612734449581,5.64960426047597,2.64048143697042,-4.68779843484001,0.896818317498234,0.38664003),
c("GO:0051302","regulation of cell division",1.05604291221993,-5.46060184840745,5.48988071528277,2.27875360095283,-2.46095888481802,0.974327411517944,0.0413688),
c("GO:0051304","chromosome separation",0.0614628149969269,6.10697750199819,-2.59381733603878,1.07918124604762,-7.02157224147822,0.771400577308842,0.49356138),
c("GO:0051783","regulation of nuclear division",0.810191652232218,0.199293452523922,-0.247925252670924,2.16435285578444,-6.7585959404082,0.884106543002424,0.04011424),
c("GO:0051988","regulation of attachment of spindle microtubules to kinetochore",0.100575515449517,-2.30052874803369,4.20585935326665,1.27875360095283,-2.48506921401123,0.756473408976317,0.47942313),
c("GO:0052695","cellular glucuronidation",0.111750572721685,-4.51931792114201,-6.0324245188599,1.32221929473392,-1.65039089804791,0.908901809807186,0.36691171),
c("GO:0071103","DNA conformation change",0.508465105883668,2.13888954564384,6.17245656725513,1.96378782734556,-5.70702639051271,0.870218117152348,0.31834856),
c("GO:0071459","protein localization to chromosome, centromeric region",0.139688215902106,-6.51111933095563,-2.86166680042405,1.41497334797082,-1.85729759822865,0.954956306196119,0.16984106),
c("GO:1902930","regulation of alcohol biosynthetic process",0.268201374532044,6.71297028853453,2.16904265265845,1.69019608002851,-3.11988504058249,0.949146572877424,0.10800743),
c("GO:2001251","negative regulation of chromosome organization",0.441414762250657,0.895073752345233,-0.455448730007491,1.90308998699194,-1.68973594059807,0.815814183463218,0.4908475));

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


chosen <- read_excel("/scratch/trcanmed/DE_RNASeq/local/share/data/chosen_for_all_revigo.xlsx", 
                     sheet = "LMO_LMH_up")
chosen$chosen <- "YES"
#chosen <- rbind(chosen , c("muscle organ development", "YES"))

one.data <- merge(one.data, chosen, by = "description", all.x = TRUE)
one.data$chosen[is.na(one.data$chosen)] <- "NO"
one.data <- one.data[, c(2,1,3,4,5,6,7,8,9,10)]
rownames(one.data) <- one.data$term_ID
one.data$new_log_size <- NA

for (i in rownames(one.data)) {
  if (one.data[i, "log_size"] < 1.96) {
    one.data[i, "new_log_size"] <- 3
  } else if (one.data[i, "log_size"] > 1.96 & one.data[i, "log_size"] < 2.5 | one.data[i, "log_size"] == 1.96) {
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
p1 <- p1 + scale_size_continuous(labels = c(0.95, 1.96, 3.0), breaks = c(3, 5, 8)) + theme_bw()
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
ggsave("revigo_lmo_lmh_up.pdf", width = 120, height = 107, useDingbats=FALSE, units = "mm")
p1
dev.off()
