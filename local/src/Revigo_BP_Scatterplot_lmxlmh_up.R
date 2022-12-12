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
revigo.data <- rbind(c("GO:0000280","nuclear division",1.81594680672738,2.87374833340027,6.31111390052016,2.51321760006794,-26.0192601045156,0.826934789766707,0.0035574),
c("GO:0000956","nuclear-transcribed mRNA catabolic process",0.575515449516679,-5.57932456826189,1.3624021747429,2.01703333929878,-7.13141994621129,0.787671601446143,0.34931451),
c("GO:0001522","pseudouridine synthesis",0.100575515449517,-7.02145839525761,3.61263564781323,1.27875360095283,-2.14978958263736,0.903399164258104,0.21531877),
c("GO:0002227","innate immune response in mucosa",0.128513158629938,-4.21704865378776,-4.42189649858951,1.38021124171161,-1.49379688648524,0.955829317683886,0.43552557),
c("GO:0006260","DNA replication",1.10633066994468,-5.77116404064579,1.92681862337755,2.29885307640971,-22.1590299565717,0.828760252899861,0.0033528),
c("GO:0006275","regulation of DNA replication",0.765491423143544,-6.96899820608278,-3.76689959633286,2.13987908640124,-5.55420552469732,0.93032149066029,0.27409895),
c("GO:0006297","nucleotide-excision repair, DNA gap filling",0.0335251718165056,-4.97606657726651,0.081327661109831,0.845098040014257,-1.49379688648524,0.834517245728671,0.47859407),
c("GO:0006325","chromatin organization",3.4083924680114,1.60057812246575,6.11846926689041,2.78604121024255,-11.9401711698425,0.910112862611965,0.22001979),
c("GO:0006353","DNA-templated transcription termination",0.145275744538191,-6.88176691839461,2.34406446488032,1.43136376415899,-1.95261586406018,0.899729568471444,0.26629341),
c("GO:0006369","termination of RNA polymerase II transcription",0.0726378722690954,-7.29486932733891,2.51781224787517,1.14612803567824,-3.0563428762315,0.904813331817324,0.20934603),
c("GO:0006614","SRP-dependent cotranslational protein targeting to membrane",0.0949879868134324,6.55355267916867,1.19227872256279,1.25527250510331,-7.98840431191851,0.931356833256381,0.00260933),
c("GO:0007007","inner mitochondrial membrane organization",0.206738559535118,3.63766394866848,5.42413100762559,1.57978359661681,-1.32307849936103,0.890435891794467,0.29023333),
c("GO:0007059","chromosome segregation",1.63714589037269,5.87863522241216,3.45943465981826,2.46834733041216,-26.8597421915997,0.766645015455795,0),
c("GO:0007077","mitotic nuclear membrane disassembly",0.0391127004525898,4.41865105669187,4.70708890560429,0.903089986991944,-2.05059310103278,0.752157306537735,0.47490697),
c("GO:0007292","female gamete generation",0.843716824048723,5.73048789439074,-1.51010429450394,2.18184358794477,-2.27451274375326,0.929605459985723,0),
c("GO:0009200","deoxyribonucleoside triphosphate metabolic process",0.0502877577247583,-3.99134226754809,4.98884241007089,1,-3.02431438553456,0.92263346475376,0.17823926),
c("GO:0009451","RNA modification",0.988992568586914,-6.23634189131701,3.13692740476547,2.25042000230889,-1.34059394738213,0.890302610073432,0.31858408),
c("GO:0016072","rRNA metabolic process",1.40805721629323,-5.76597446137749,3.00581755902588,2.40312052117582,-7.05517038561118,0.863407547967649,0.2801338),
c("GO:0019080","viral gene expression",0.290551489076382,-2.3899704154271,8.6386529921407,1.72427586960079,-7.08263600799384,0.996005962227381,0),
c("GO:0022613","ribonucleoprotein complex biogenesis",2.49203777169358,1.88725838475074,7.26769548901543,2.65030752313194,-9.10856409735512,0.894986482678361,0.20104668),
c("GO:0031055","chromatin remodeling at centromere",0.0447002290886741,0.762543654501582,4.1907100545088,0.954242509439325,-13.7683477289923,0.924531297002653,0.13866768),
c("GO:0031123","RNA 3'-end processing",0.6034530926971,-6.0485555724116,3.95117342252863,2.03742649794062,-1.30301696301539,0.87481193720032,0.47885017),
c("GO:0031640","killing of cells of another organism",0.206738559535118,-4.10405525048153,-5.35689161689267,1.57978359661681,-1.30391371277847,0.984966310408338,0.37578685),
c("GO:0032886","regulation of microtubule-based process",1.40246968765715,-3.92995555083807,8.35417546128689,2.40140054078154,-1.9892794957791,0.973571269699561,0.04280184),
c("GO:0042254","ribosome biogenesis",1.68184611946136,1.56024444837079,7.44267114546954,2.48000694295715,-8.95129810954029,0.871812612018782,0.41746946),
c("GO:0042276","error-prone translesion synthesis",0.0558752863608426,-5.56279911709113,0.133361571775304,1.04139268515823,-2.53673727129018,0.821111475252516,0.49897757),
c("GO:0042398","cellular modified amino acid biosynthetic process",0.229088674079455,-8.22160047152985,-1.05313148701401,1.6232492903979,-1.46788954894494,0.941084568407914,0.33803552),
c("GO:0042770","signal transduction in response to DNA damage",0.771078951779628,-3.20589341997475,-2.44077127225973,2.14301480025409,-3.9530696620242,0.873167999455989,0.44091118),
c("GO:0042775","mitochondrial ATP synthesis coupled electron transport",0.519640163155836,-1.1377644721779,8.74945766549202,1.9731278535997,-2.42001704420421,0.92622509712078,0.09297333),
c("GO:0043628","regulatory ncRNA 3'-end processing",0.106163044085601,-6.12021373033076,4.38727121810655,1.30102999566398,-1.54205114612257,0.889082396688405,0.21634428),
c("GO:0044364","disruption of cells of another organism",0.041906464770632,-3.68789925161383,-5.48412779610939,0.903089986991944,-1.30391371277847,0.989329585206785,0.33081295),
c("GO:0045023","G0 to G1 transition",0.0167625859082528,5.85582714828657,3.19292888457488,0.602059991327962,-2.32660629123618,0.83575950769811,0.44348514),
c("GO:0045605","negative regulation of epidermal cell differentiation",0.0726378722690954,-0.752788854872795,-5.88830574889646,1.14612803567824,-2.81105788141478,0.928019089172925,0.0314378),
c("GO:0045682","regulation of epidermis development",0.363189361345477,-1.10295837259435,-6.55156193957715,1.81954393554187,-1.31584262909079,0.974438172911217,0.26593307),
c("GO:0045787","positive regulation of cell cycle",1.97798513717383,1.49827142252472,-5.91344296751248,2.55022835305509,-12.3756269360325,0.743313012744154,0),
c("GO:0046471","phosphatidylglycerol metabolic process",0.184388444990781,-2.74829891378654,7.09225391075936,1.53147891704226,-1.46845894723707,0.956422181231688,0.40296563),
c("GO:0048285","organelle fission",1.95563502262949,2.89568183545397,5.90565179926122,2.54530711646582,-24.1579039919154,0.871631178793252,0.37233798),
c("GO:0048665","neuron fate specification",0.178800916354696,6.35687021353132,-0.611404600221912,1.51851393987789,-1.60884480349048,0.984671529277251,0.12040871),
c("GO:0050000","chromosome localization",0.458177348158909,5.64997827269408,1.13109079907317,1.91907809237607,-5.80149195436448,0.949976127646982,0.16453509),
c("GO:0050829","defense response to Gram-negative bacterium",0.474939934067162,-4.11656606823953,-4.11626407294331,1.93449845124357,-1.94960758958804,0.951124693747117,0.48771656),
c("GO:0051052","regulation of DNA metabolic process",2.99491534894116,-7.21375524849848,-3.31815651666331,2.72997428569956,-6.56578761605521,0.949143685642251,0.04717156),
c("GO:0051276","chromosome organization",2.43616248533274,2.6516932579456,5.68568117720859,2.64048143697042,-14.7725282021257,0.868799887959448,0.38664003),
c("GO:0051302","regulation of cell division",1.05604291221993,0.748609833897075,0.784936430788898,2.27875360095283,-3.88019204009741,0.974308701247368,0.0413688),
c("GO:0051438","regulation of ubiquitin-protein transferase activity",0.324076660892887,3.4944147635077,-4.70991262108022,1.77085201164214,-1.35289611543266,0.975402541857666,0.46106067),
c("GO:0051726","regulation of cell cycle",6.27479465832262,2.76902451828912,0.739019707379006,3.05076631123304,-3.37107189446829,0.968858738353463,0.05939145),
c("GO:0051783","regulation of nuclear division",0.810191652232218,0.982250490239902,-3.31278771648744,2.16435285578444,-11.5482298654284,0.903566854029217,0.04011424),
c("GO:0060249","anatomical structure homeostasis",1.25160641448287,4.69880867091007,-0.72635513413836,2.35218251811136,-1.37604960434702,0.98964406165366,0.1459994),
c("GO:0061844","antimicrobial humoral immune response mediated by antimicrobial peptide",0.385539475889814,-3.86653914947355,-4.25306618455059,1.84509804001426,-2.39387955441391,0.950330657326216,0.07823026),
c("GO:0070125","mitochondrial translational elongation",0.0279376431804213,-7.35935128807138,0.803123577941111,0.778151250383644,-5.58365016575794,0.889604461063442,0.26601885),
c("GO:0070317","negative regulation of G0 to G1 transition",0.0335251718165056,0.915073034618103,-5.74670893648496,0.845098040014257,-2.65697686708129,0.75270394466524,0.43525067),
c("GO:0071103","DNA conformation change",0.508465105883668,3.45925901936692,6.21021680822973,1.96378782734556,-20.9124192738021,0.839215845177389,0.31834856),
c("GO:0071824","protein-DNA complex subunit organization",1.29071911493546,0.95234510448968,6.66467130866345,2.3654879848909,-8.41138155179415,0.89448547094119,0.19445929),
c("GO:0072331","signal transduction by p53 class mediator",0.497290048611499,-0.902443212621817,-1.49454531953731,1.95424250943932,-4.33455359593377,0.918165252775678,0.03799221),
c("GO:0090305","nucleic acid phosphodiester bond hydrolysis",1.525395317651,-4.95443448315304,2.83321419059966,2.43775056282039,-2.27755600401472,0.892364252205076,0.29078755),
c("GO:0090502","RNA phosphodiester bond hydrolysis, endonucleolytic",0.519640163155836,-5.20245643002346,3.39951069805577,1.9731278535997,-2.91969026896702,0.883555322141271,0.29705407),
c("GO:0140053","mitochondrial gene expression",0.787841537687881,-6.16704697343599,6.24923311306772,2.15228834438306,-2.24019390180073,0.935587051520974,0.11727514),
c("GO:1904666","regulation of ubiquitin protein ligase activity",0.150863273174275,3.38028880647766,-4.93679789332311,1.44715803134222,-2.29313493173082,0.976790987705713,0.02172857),
c("GO:1904874","positive regulation of telomerase RNA localization to Cajal body",0.0838129295412639,2.40594171089662,-1.61992396513179,1.20411998265592,-1.60884480349048,0.973527288211116,0.12342507),
c("GO:2001020","regulation of response to DNA damage stimulus",1.74889646309437,5.08016786396543,-3.25026982296583,2.49692964807321,-1.30764536049963,0.950266476567345,0.32177574),
c("GO:2001038","regulation of cellular response to drug",4.35561292740989,4.7043679673021,-3.15393669630008,2.89486965674525,-1.32307849936103,0.957139413261742,0.05277783),
c("GO:2001251","negative regulation of chromosome organization",0.441414762250657,0.338775790652921,-3.8742814387123,1.90308998699194,-4.88823565753134,0.845225252036609,0.4908475));

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
                     sheet = "LMX_LMH_up")
chosen$chosen <- "YES"
#chosen <- rbind(chosen , c("muscle organ development", "YES"))

one.data <- merge(one.data, chosen, by = "description", all.x = TRUE)
one.data$chosen[is.na(one.data$chosen)] <- "NO"
one.data <- one.data[, c(2,1,3,4,5,6,7,8,9,10)]
rownames(one.data) <- one.data$term_ID
one.data$new_log_size <- NA

for (i in rownames(one.data)) {
  if (one.data[i, "log_size"] < 1.95) {
    one.data[i, "new_log_size"] <- 0.60
  } else if (one.data[i, "log_size"] > 1.95 & one.data[i, "log_size"] < 2.5 | one.data[i, "log_size"] == 1.95) {
    one.data[i, "new_log_size"] <- 1.95
  } else {
    one.data[i, "new_log_size"] <- 3.0
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
p1 <- p1 + scale_size_continuous(breaks = c(0.60, 1.95, 3.0)) + theme_bw()
ex <- one.data [ one.data$chosen == "YES", ];
#ex <- one.data [ one.data$dispensability < 0.15, ]
p1 <- p1 + geom_label_repel(data = ex, aes(plot_X, plot_Y, label = description), color = 'black', size = 3,
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

ggsave("revigo_lmx_lmh_up.pdf", width = 11.69, height = 8.27)
p1 <- p1 + unmute_theme
dev.off()
