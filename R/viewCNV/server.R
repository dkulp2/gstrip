#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(plyr)
library(ggplot2)
library(stats)
library(foreach)
library(zoo)
library(gridExtra)
library(tidyr)

# input <-
#   structure(list(candidate = "DEL_P0563_29", seg.chr = 20L, seg.end = 36477348L, 
#                  seg.start = 36474476L), .Names = c("candidate", "seg.chr", 
#                                                     "seg.end", "seg.start"))

shell <- "C:\\cygwin64\\bin\\bash.exe -lc"

input.dump.fn <- "d:/mccarroll/tmp/inputdump.R"

pad <- 5000
extent.offset <- 12*100/2 # FIXME: outputs should set map positions centered in their windows
#extent.offset <- 0
profile.fn <- "d:/mccarroll/gpc_wave2_batch1/profile_seq_20_100.dat.gz"
profile.fn <- '/home/dkulp/mccarroll/gpc_wave2_batch1/profile_seq_20_100.dat.gz'

data.dir <- "d:/mccarroll/cnv_seg.12.500"
irs.fn <- paste0(data.dir,"/cnv_segs.irs")
gs_dels.fn <- paste0(data.dir,"/gs_dels.txt")
probe.fn <- paste0(data.dir,"/probes.txt")

cn.segs.merged.fn <- paste0(data.dir,"/sites_cnv_segs.txt.csm.Rdata")
cnv.geno.fn <- paste0(data.dir,"/sites_cnv_segs.txt.cnvgeno.srt.gz")

log.fn <- paste0(data.dir,"/notes.Rdata")
if (file.exists(log.fn)) {
  load(log.fn)
  cat(sprintf("Loaded %s notes from %s",nrow(notes),log.fn))
} else {
  notes <- data.frame()
}

frags.header <- unlist(strsplit(system(sprintf("%s 'zcat %s | head -1", shell, profile.fn), intern=T),"\t"))
genoWin.header <- unlist(strsplit(system(sprintf("%s 'zcat %s | head -1", shell, cnv.geno.fn), intern=T),"\t"))

# load simple CNV fragments as a starting point. These are extents with a consistent call in adjacent windows.
load(cn.segs.merged.fn)
cn.segs.merged$evidence <- 'profile extents'
#cn.segs.merged$start.map <- cn.segs.merged$start.map + extent.offset
cn.segs.merged$end.map <- cn.segs.merged$end.map + extent.offset

top.deletions <- (cn.segs.merged[cn.segs.merged$cn==1,] %>% arrange(-len) %>% head(400))[['seg']]
cat(sprintf("Loaded %s CNVs from %s\n",nrow(cn.segs.merged), cn.segs.merged.fn))

irs.orig <- read.table(irs.fn,header=T,sep="\t", as.is=T)
colnames(irs.orig) <- c('seg','chr','start','end','Pval','nprobes','# Samples','Lower Pval','Lower # Samples','Higher Pval', 'Higher # Samples')
cat(sprintf("Loaded %s IRS rows from %s\n", nrow(irs.orig),irs.fn))

gs.dels.orig <- read.table(gs_dels.fn, header=TRUE, sep="\t", stringsAsFactors = FALSE)
colnames(gs.dels.orig) <- c('.id','seg','chr','start.map','end.map','cn','cq','paired.reads')
gs.dels.orig$evidence <- 'gstrip dels'
gs.dels.orig$copy.number <- addNA(as.factor(gs.dels.orig$cn))
top.gs.deletions <- (arrange(subset(gs.dels.orig, end.map-start.map < 8000), -(end.map-start.map))[['seg']] %>% unique)
cat(sprintf("Loaded %s GStrip deletions from %s\n", nrow(gs.dels.orig),gs_dels.fn))

probes.orig <- read.table(probe.fn)
colnames(probes.orig) <- c('seg','chr','start.map','end.map')
cat(sprintf("Loaded %s probes from %s\n",nrow(probes.orig),probe.fn))

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  # load the selection options for candidate sites
  updateSelectInput(session, 'candidate', choices=c("Choose Known Deletion"="",top.gs.deletions))
  
  # load the selection options for saved sites
  updateSelectInput(session, 'saved.site', choices=c("Choose Preset"='',notes$name))

  # load the values for all sample names
  updateSelectizeInput(session, 'seg.sample', choices=unique(cn.segs.merged$.id))

  # if the user selects a saved site, load those coords
  observe({
    input$reset
    
    # save the latest selected site for interactive debug
    input.dump <- reactiveValuesToList(input)
    dump('input.dump',file=input.dump.fn)
    
    if (input$saved.site != "") {
      note <- notes[notes$name == input$saved.site,]
      updateTextInput(session,'seg.chr', value=note$chr)
      updateNumericInput(session,'seg.start', value=note$start.map)
      updateNumericInput(session,'seg.end', value=note$end.map)
      updateSelectizeInput(session, 'seg.sample', selected=unlist(strsplit(note$sample,',',fixed=TRUE)))

      # must call special added hook to set value of textarea because there is no builtin R func
      session$sendCustomMessage(type = 'setnote', message = list(note = note$note))
    }
  })
  
  # prev saved
  observe({
    input$saved.prev
    idx <- which(notes$name == isolate(input$saved.site))
    if (length(idx)>0 && idx > 1) { updateSelectInput(session,'saved.site',selected=notes$name[idx-1]) }
  })
  
  # next saved
  observe({
    input$saved.next
    idx <- which(notes$name == isolate(input$saved.site))
    if (length(idx)>0 && idx < nrow(notes)) { updateSelectInput(session,'saved.site',selected=notes$name[idx+1]) }
  })

    # if the candidate site changes, then set the samples, chromosome, start and end positions
  observe({
    if (input$candidate != "") {
      gs.row <- gs.dels.orig[gs.dels.orig$seg == input$candidate,]
      seg.sample <- gs.row[,'.id']
      seg.chr <- gs.row[1,'chr']
      seg.start <- gs.row[1,'start.map']
      seg.end <- gs.row[1,'end.map']
      updateTextInput(session,'seg.chr', value=seg.chr)
      updateNumericInput(session,'seg.start', value=seg.start)
      updateNumericInput(session,'seg.end', value=seg.end)
      updateSelectizeInput(session, 'seg.sample', selected=seg.sample)
    }
  })
  
  # prev candidate
  observe({
    input$candidate.prev
    idx <- which(top.gs.deletions == isolate(input$candidate))
    if (length(idx)>0 && idx > 1) { updateSelectInput(session,'candidate',selected=top.gs.deletions[idx-1]) }
  })
  
  # next candidate
  observe({
    input$candidate.next
    idx <- which(top.gs.deletions == isolate(input$candidate))
    if (length(idx)>0 && idx < length(top.gs.deletions)) { updateSelectInput(session,'candidate',selected=top.gs.deletions[idx+1]) }
  })
  
  # zoom out 1.5x
  observe({
    input$zoom.out
    cat("zoom out\n")
    incr <- isolate(floor(.25 * (input$seg.end-input$seg.start)))
    updateNumericInput(session,'seg.start', value=max(0,isolate(input$seg.start - incr )))
    updateNumericInput(session,'seg.end', value=isolate(input$seg.end + incr))
  })
  
  observe({
    input$delete
    isolate({
      notes <<- notes[notes$name != input$saved.site,]
      save(notes, file=log.fn)
      updateSelectInput(session, 'saved.site', choices=c("Choose Preset"='',notes$name))
    })
  })
  
  # prepend site with notes to log file.
  observe({
    input$save.button
    isolate({
      if (input$note != "" && length(input$seg.sample) > 0) {
        cat("save\n")
        note.name <- sprintf("chr%s %sMb (%skb) %s",
                             input$seg.chr, format(input$seg.start/1000000,digits=4,nsmall=2),
                             format((input$seg.end-input$seg.start)/1000,digits=4,nsmall=2),
                             regmatches(input$note, regexpr("[^\n]{1,50}",input$note,perl=TRUE)))
        new.note <- data.frame(name=note.name, chr=input$seg.chr, start.map=input$seg.start, end.map=input$seg.end,
                               sample=paste(input$seg.sample,collapse=','),note=input$note,
                               stringsAsFactors = FALSE)
        if (nrow(notes)==0) { notes <<- new.note } else { 
          if (any(notes$name==note.name)) {
            notes <<- rbind(new.note, notes[-which(notes$name==note.name),])
          } else {
            notes <<- rbind(new.note, notes)  
          }
        }
        notes <<- unique(notes)
        save(notes, file=log.fn)
        updateSelectInput(session, 'saved.site', choices=c("Choose Preset"='',notes$name), selected=note.name)
      }
    })
  })

  # a subset of the extents
  cn.segs <- reactive({
    x <- cn.segs.merged[cn.segs.merged$chr == input$seg.chr & 
                          cn.segs.merged$end.map > input$seg.start-input$pad & 
                          cn.segs.merged$start.map < input$seg.end+input$pad & 
                          cn.segs.merged$cn != 2 & 
                          cn.segs.merged$len > input$min.cnv.len & 
                          !is.na(cn.segs.merged$cn),]
    x$target <- x$.id %in% input$seg.sample
    return(x)
  })
  
  irs <- reactive({
    subset(irs.orig, chr==input$seg.chr & end > input$seg.start-input$pad & start < input$seg.end+input$pad & nprobes > 0)
  })
  
  probes <- reactive({
    subset(probes.orig, chr==input$seg.chr & end.map > input$seg.start-input$pad & start.map < input$seg.end+input$pad)
  })
  
  gs.dels <- reactive({
    x <- subset(gs.dels.orig, chr==input$seg.chr & end.map > input$seg.start-input$pad & start.map < input$seg.end+input$pad & paired.reads >= input$paired.reads)
    x$target <- x$.id %in% input$seg.sample
    return(x)  
  })
  
  # format gs.dels, cn.segs, probes, and irs into a single data.frame for ggplot
  cnv.disp <- reactive({
    gs.dels <- gs.dels()
    cn.segs <- cn.segs()
    irs <- irs()
    probes <- probes()
    
    disp.cols <- c('.id','cn','copy.number','chr','start.map','end.map','seg', 'target', 'evidence')
    gs.del.sel <- gs.dels[,c(disp.cols,'paired.reads')] 
    cn.segs.sel <- cn.segs[,disp.cols]
    if (nrow(cn.segs.sel)>0) { cn.segs.sel$paired.reads <- NA }
    print(cn.segs.sel)
    disp <- rbind(gs.del.sel, cn.segs.sel)
    
    disp2 <- rbind(disp, 
                   data.frame(.id='Probe',cn=NA,copy.number=NA,chr=probes$chr,start.map=probes$start.map,end.map=probes$end.map,seg=probes$seg,target=FALSE,evidence="Array",paired.reads=NA))
    
    disp3 <- join(disp2, subset(irs, select=-c(chr,start,end)), by="seg", type="left")
    disp3$Pval <- ifelse(is.na(disp3$Pval),1,disp3$Pval)
    disp3$Phred <- ifelse(disp3$Pval==0, 100, -log10(disp3$Pval)*10)
    return(disp3)
  })

  winGeno <- reactive({
    # load windowed genotypes
    wg.fn <- basename(tempfile("wg"))
    tabix <- sprintf("tabix -h %s %s:%s-%s > /tmp/%s", cnv.geno.fn, input$seg.chr, input$seg.start-input$pad, input$seg.end+input$pad, wg.fn)
    system(sprintf("%s '%s'", shell, tabix))
    wg <- read.table(sprintf("C:\\cygwin64\\tmp\\%s",wg.fn), header=FALSE, sep="\t", as.is=TRUE, check.names=FALSE)
    colnames(wg) <- genoWin.header
    wg$cn <- as.factor(wg$cn)
    wg$pos <- wg$start + (wg$end - wg$start)/2
    return(wg)
  })
  
  frags <- reactive({
    # load fragment data
    frag.fn <- basename(tempfile("cnv"))
    tabix <- sprintf("tabix -h %s %s:%s-%s > /tmp/%s", profile.fn, input$seg.chr, input$seg.start-input$pad, input$seg.end+input$pad, frag.fn)
    system(sprintf("%s '%s'", shell, tabix))
    frags <- read.table(sprintf("C:\\cygwin64\\tmp\\%s",frag.fn), header=FALSE, sep="\t", as.is=TRUE, check.names=FALSE)
    colnames(frags) <- frags.header
    return(frags)
  })
  
  # returns a data.frame for plotting that is a row per sample and expected counts
  expected <- reactive({
    frags <- frags()
    sample.count <- ncol(frags) - 6
    exp <- foreach(i=1:sample.count) %do% {
      obs_exp <- strsplit(frags[,6+i], ',')
      as.numeric(unlist(lapply(obs_exp, function(x) { x[2] })))
    }
    names(exp) <- colnames(frags[1,7:(6+sample.count)])
    exp.df <- as.data.frame(exp)
    colnames(exp.df) <- names(exp)
    exp.df$pos <- as.integer(1:nrow(exp.df))
    exp2 <- gather(exp.df, sample, exp, -pos)
    exp2$start.map <- frags[exp2$pos,'START']
    exp2$sample.type <- factor(ifelse(exp2$sample %in% input$seg.sample, 'target','other'))
    return(exp2)
  })
  
  # observed/expected ratio
  oer <- reactive({
    frags <- frags()
    sample.count <- ncol(frags) - 6
    # returns a list of ratios for each window
    obs_exp_ratio <- foreach(i=1:sample.count) %do% {
      # each column is comma-separated observed,expected
      obs_exp <- strsplit(frags[,6+i], ',')
      obs <- as.integer(unlist(lapply(obs_exp, function(x) { x[1] })))
      exp <- as.numeric(unlist(lapply(obs_exp, function(x) { x[2] })))
      
      # window obs and exp
      obs_mean <- rollapply(obs, width=input$win.size, sum)
      exp_mean <- rollapply(exp, width=input$win.size, sum)
      obs_mean / exp_mean
    }
    
    names(obs_exp_ratio) <- colnames(frags[1,7:(6+sample.count)])
    obs_exp_ratio.df <- as.data.frame(obs_exp_ratio)
    colnames(obs_exp_ratio.df) <- names(obs_exp_ratio)

    obs_exp_ratio.df$site.median = apply(obs_exp_ratio.df, 1, median, na.rm=TRUE)
    obs_exp_ratio.df$pos <- as.integer(1:nrow(obs_exp_ratio.df))
    
    oer <- gather(obs_exp_ratio.df, sample, ratio, -pos)
    oer$start.map <- frags[oer$pos,'START'] + (frags[oer$pos + input$win.size - 1,'END'] - frags[oer$pos,'START'])/2
    oer$sample.type <- ifelse(oer$sample == 'site.median','median', ifelse(oer$sample %in% input$seg.sample, 'target','other'))
    return(oer)  
  })
  
  xbounds <- reactive({
    xlim(min(oer()$start.map,cnv.disp()$start.map),max(oer()$start.map,cnv.disp()$end.map))
  })
  
  cnvPlot <- reactive({
    cnv.disp <- cnv.disp()
    if (input$show_readPairs) {
      plt <- ggplot(cnv.disp, aes(x=start.map, xend=end.map, y=.id, yend=.id, color=paired.reads, size=target)) + scale_color_gradient2() 
    } else {
      plt <- ggplot(cnv.disp, aes(x=start.map, xend=end.map, y=.id, yend=.id, color=Phred, size=target)) + scale_color_gradient(limits=c(0,100))

    }
    return(plt + geom_segment() + facet_grid(evidence ~ ., scales="free", space="free") + geom_point(data=subset(cnv.disp, evidence=='Array'), size=5, color='black') + guides(size="none") + theme(axis.title.x = element_blank()) + xbounds())
  })
  
  winGenoPlot <- reactive({
    df <- winGeno()
    ggplot(subset(df, sample %in% intersect(gs.dels()$.id,input$seg.sample)), aes(x=pos, group=sample, color=cnq, y=cn)) + scale_color_gradient(limits=c(0,100)) + geom_point(size=5, shape=15) + facet_grid(sample ~ ., scale="free") + xbounds() + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,as.numeric(input$winGeno_nudge_L)),'in'))

  })
  
  fragPlot <- reactive({
    df <- oer()
    ggplot(df, aes(x=start.map, y=ratio, group=sample)) + geom_step(color='grey') +
     geom_step(data=df[df$sample.type=='median',], color='black') + 
      geom_step(data=df[df$sample.type=='target',], color='red') +
      theme(axis.title.x = element_blank(), plot.margin=unit(c(0,1,0,0.5),'inches')) +
      xbounds()
  })
  
  output$expPlot <- renderPlot({
    expected <- expected()
    colors <- c(target='#FF0000',other='#222222')
    alpha <- c(target=1,other=0.2)
    ggplot(expected, aes(x=start.map,y=exp,group=sample,color=sample.type,alpha=sample.type)) + geom_line() +
      theme(axis.title.x = element_blank(), plot.margin=unit(c(0,1,0,0.5),'inches')) + 
      scale_color_manual(name = "sample.type",values = colors) +
      scale_alpha_manual(name = "sample.type",values = alpha) +
      xbounds() + guides(color=FALSE,alpha=FALSE)
  })
  
  output$seg.size <- renderText({ paste0(format((input$seg.end-input$seg.start)/1000,digits=4,nsmall=2),"kb") })
  
  output$allThree <- renderPlot({
    plots <- list()
    if (input$show_cnv) { plots$cnv = cnvPlot() }
    if (input$show_winGeno) { plots$winGeno = winGenoPlot() }
    if (input$show_frag) { plots$frag = fragPlot() }
    grid.arrange(do.call(arrangeGrob,c(plots,list(ncol=1))))
  })
})
