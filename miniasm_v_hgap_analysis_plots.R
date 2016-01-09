#!/usr/bin/env Rscript

library(ggplot2)

a=read.csv(file="miniasm_v_hgap.assembly_stats.tsv", sep="\t", header=TRUE)

t = theme(text=element_text(size=15), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1,size=10))

ggplot(data=a) +
  geom_bar(aes(x=sample,y=total_length,fill=assembler), stat="identity", width=0.75, position=position_dodge()) +
  t +
  xlab("Sample") +
  ylab("Total bases in assembly")
ggsave("miniasm_v_hgap_analysis.plot.total_bases_compare.png")



ggplot(data=a) +
  geom_bar(aes(x=sample,y=number,fill=assembler), stat="identity", width=0.75, position=position_dodge()) +
  t +
  xlab("Sample") +
  ylab("Number of contigs")
ggsave("miniasm_v_hgap_analysis.plot.number_of_contigs.png")



a=read.csv(file="miniasm_v_hgap.hit_length_compare.tsv", sep="\t", header=TRUE)

ggplot(data=a) +
  geom_boxplot(aes(x=Sample,y=Length_ratio, fill=Assembler)) +
  t +
  xlab("Sample") +
  ylab("Assembly hit length (% of reference hit length)")
ggsave("miniasm_v_hgap_analysis.plot.hit_length_ratio.png")




library(reshape)
a=read.csv(file="miniasm.resources.tsv", sep="\t", header=TRUE)
a=melt(a, id.vars="Sample")
time_data = subset(a, variable %in% c("cpu_time", "wall_clock_time"))
colnames(time_data)[2] = "Type"
colnames(time_data)[3] = "Time"
time_data$Time = time_data$Time / 60
ggplot(data=time_data) +
  geom_bar(aes(x=Sample, y=Time, fill=Type), stat="identity", width=0.75, position=position_dodge()) +
  t +
  xlab("Sample") +
  ylab("Time (minutes)")
ggsave("miniasm.resources.run_time.png")


mem = subset(a, variable %in% c("max_memory"))
ggplot(data=mem) +
  geom_bar(aes(x=Sample, y=value), stat="identity", width=0.75) +
  t +
  xlab("Sample") +
  ylab("Peak memory (GB)")
ggsave("miniasm.resources.memory.png")

