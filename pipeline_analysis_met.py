#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import math
import pandas as pd



    
def Diversity_methylation_program(ArgsVal):
    global args
    description = """ Study of methylation diversity """
    parser = argparse.ArgumentParser(prog="methylation diversity",description=description)
    parser.add_argument("-P", "--parameter", default="parameter_PARENT.txt")
    parser.add_argument("-C", "--clone", default="parameter_CLONE.txt")
    parser.add_argument("-D", "--directory", default="/")
    parser.add_argument("-V", "--vcftools", default="vcftools")
    parser.add_argument("-N", "--Pop_nb", default=7)
    parser.add_argument("-b", "--bind", default="T")
    parser.add_argument("-cv", "--convert", default="T")
    parser.add_argument("-k", "--kinship", default="T")
    parser.add_argument("-p", "--population_analysis", default="T")
    parser.add_argument("-cs", "--conservation", default="T")
    parser.add_argument("-n", "--name_analysis", default="CG")
    parser.add_argument("-T", "--met_treshold", default=50)
    parser.add_argument("-dl", "--dl_windows_size", default=30000)

    

    args = parser.parse_args(ArgsVal)
    parameter=str(args.parameter)
    nb_pop=int(args.Pop_nb)
    met_treshold=int(args.met_treshold)
    dl_windows_size=int(args.dl_windows_size)
    clone=str(args.clone)
    pop=str(args.name_analysis)
    directory=str(args.directory)
    vcftools=str(args.vcftools)
    bind=str(args.bind)
    convert=str(args.convert)
    population_analysis=str(args.population_analysis)
    kinship=str(args.kinship)
    conservation=str(args.conservation)
    Parameter=open(str(parameter),"r")
    Parameter=Parameter.read()
    Parameter=Parameter.split('#')
    Clone=open(str(clone),"r")
    Clone=Clone.read()
    Clone=Clone.split('#')
    

    # 1 bind of all bedgrap files
    if bind=="T":
        bedgraph=Parameter[1]
        bedgraph = bedgraph.split("\n")
        bedgraph=bedgraph[1:]
        pos_parent=[]
        id_ind=[]
        for file in bedgraph:
            if file!="":
                id_ind.append(file)
                name=str(file)+".bedGraph"
                bedgraph2 = open(str(name), "r")
                bedgraph2=bedgraph2.read()
                bedgraph2 = bedgraph2.split("\n")
                exec(str(file)+"_dict={}")
                for line in bedgraph2[1:]:
                    if line!="":
                        line_splitted = line.split("\t")
                        pos=line_splitted[0]+"_"+str(line_splitted[1])
                        if pos not in pos_parent:pos_parent.append(pos)
                        var="0/0"
                        if int(line_splitted[3])<met_treshold:var="1/1"
                        exec(str(file)+"_dict[str(pos)]=var")
        print("number of called sites:"+str(len(pos_parent)))
        pos_parent.sort()


        fileout=open("Genotype_parent_"+str(pop)+".csv","w")
        written_line="#rs\tallele\tchrom\tpos"
        for ind in id_ind:
            written_line=str(written_line)+"\t"+str(ind)
        fileout.write(written_line)
        nb_fixed=0
        nb_fixed_met=0
        nb_fixed_unmet=0
        for pos in pos_parent:
            name=str(pos)
            pos=pos.split("_")
            written_line="\n"+str(name)+"\tA/T\t"+str(pos[0])+"\t"+str(pos[1])
            pos2=[]
            for ind in id_ind:
                var="NA"
                exec("if name in "+str(ind)+"_dict.keys():var=str("+str(ind)+"_dict[str(name)])")
                if var!="NA":
                    pos2.append(var)
                    if var=="0/0": var="A"
                    if var=="1/1": var="T"
                written_line=str(written_line)+"\t"+str(var)
            fileout.write(written_line)
            pos2=list(set(pos2))
            if len(pos2)==1:
                nb_fixed=nb_fixed+1
                if pos2[0]=="0/0":nb_fixed_met=nb_fixed_met+1
                if pos2[0]=="1/1":nb_fixed_unmet=nb_fixed_unmet+1
        print("number of fixed sites:"+str(nb_fixed))
        print("number of fixed methylated sites:"+str(nb_fixed_met))
        print("number of fixed unmethylated sites:"+str(nb_fixed_met))
        fileout.close()

        bedgraph_clone= Clone[1]
        bedgraph_clone = bedgraph_clone.split("\n")
        bedgraph_clone=bedgraph_clone[1:]
        pos_clone=[]
        pos_all=[]
        id_ind_clone=[]
        for file in bedgraph_clone:
            if file!="":
                id_ind.append(file)
                name=str(file)+".bedGraph"
                bedgraph2 = open(str(name), "r")
                bedgraph2=bedgraph2.read()
                bedgraph2 = bedgraph2.split("\n")
                exec(str(file)+"_dict_clone={}")
                for line in bedgraph2[1:]:
                    if line!="":
                        line_splitted = line.split("\t")
                        pos=line_splitted[0]+"_"+str(line_splitted[1])
                        if pos not in pos_clone:pos_clone.append(pos)
                        if pos in pos_parent:pos_all.append(pos)
                        var="0/0"
                        if int(line_splitted[3])<met_treshold:var="1/1"
                        exec(str(file)+"_dict_clone[str(pos)]=var")
        print("number of called sites (clone):"+str(len(pos_clone)))
        pos_clone.sort()
        


        fileout=open("Genotype_clone_"+str(pop)+".csv","w")
        written_line="#rs\tallele\tchrom\tpos"
        for ind in id_ind:
            written_line=str(written_line)+"\t"+str(ind)
        fileout.write(written_line)
        nb_fixed=0
        nb_fixed_met=0
        nb_fixed_unmet=0
        for pos in pos_clone:
            name=str(pos)
            pos=pos.split("_")
            written_line="\n"+str(name)+"\tA/T\t"+str(pos[0])+"\t"+str(pos[1])
            pos2=[]
            for ind in id_ind:
                var="NA"
                exec("if name in "+str(ind)+"_dict_clone.keys():var=str("+str(ind)+"_dict_clone[str(name)])")
                if var!="NA":
                    pos2.append(var)
                    if var=="0/0": var="A"
                    if var=="1/1": var="T"
            fileout.write(written_line)
            pos2=list(set(pos2))
            if len(pos2)==1:
                nb_fixed=nb_fixed+1
                if pos2[0]=="0/0":nb_fixed_met=nb_fixed_met+1
                if pos2[0]=="1/1":nb_fixed_unmet=nb_fixed_unmet+1
        print("number of fixed sites (clone):"+str(nb_fixed))
        print("number of fixed methylated sites (clone):"+str(nb_fixed_met))
        print("number of fixed unmethylated sites (clone):"+str(nb_fixed_met))
        fileout.close()

    # 2 convert to vcf
    if convert=="T":
        x2=str(pop)+"_parent.vcf" 
        fichier=open(str(x2), "a")
        fichier.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for ind in id_ind:
            fichier.write("\t"+str(ind))
        annotation=open("Genotype_parent_"+str(pop)+".csv", "r")
        annotation=annotation.read()
        annotation = annotation.split("\n")
        for line in annotation[1:]:
            line_splitted = line.split("\t")
            chrom=line_splitted[2]
            pos=line_splitted[3]
            alt=line_splitted[1]
            alt=alt.split("/")
            ref=alt[0]
            alt=alt[1]
            fichier.write("\n"+str(chrom)+"\t"+str(pos)+"\t.\t"+str(ref)+"\t"+str(alt)+"\t1000\t.\t.\tGT")
            for i in line_splitted[4:]:
                if str(i)==str(ref):fichier.write("\t0/0")
                elif str(i)==str(alt):fichier.write("\t1/1")
                else : fichier.write("\t./.")
    
        fichier.close()
        x2=str(pop)+"_clone.vcf" 
        fichier=open(str(x2), "a")
        fichier.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for ind in id_ind:
            fichier.write("\t"+str(ind))
        annotation=open("Genotype_clone_"+str(pop)+".csv", "r")
        annotation=annotation.read()
        annotation = annotation.split("\n")
        for line in annotation[1:]:
            line_splitted = line.split("\t")
            chrom=line_splitted[2]
            pos=line_splitted[3]
            alt=line_splitted[1]
            alt=alt.split("/")
            ref=alt[0]
            alt=alt[1]
            fichier.write("\n"+str(chrom)+"\t"+str(pos)+"\t.\t"+str(ref)+"\t"+str(alt)+"\t1000\t.\t.\tGT")
            for i in line_splitted[4:]:
                if str(i)==str(ref):fichier.write("\t0/0")
                elif str(i)==str(alt):fichier.write("\t1/1")
                else : fichier.write("\t./.")
    
        fichier.close()
        x2=str(pop)+"_parent_sorted.vcf" 
        fichier=open(str(x2), "a")
        fichier.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for ind in id_ind:
            fichier.write("\t"+str(ind))
        file_sync = pd.read_csv(str(pop)+"_parent.vcf" ,sep="\t")
        file_sync=file_sync.sort_values(by=["#CHROM","POS"])
        for i in range(0,len(file_sync.index)):
            seq=file_sync.iloc[i]
            x=str(seq[0])
            for j in range(1,9):x=str(x)+"\t"+str(seq[j])
            for j in range(9+4*pop,9+4*pop+4):x=str(x)+"\t"+str(seq[j])
            fichier.write("\n"+str(x))
        fichier.close()
        x2=str(pop)+"_clone_sorted.vcf" 
        fichier=open(str(x2), "a")
        fichier.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for ind in id_ind:
            fichier.write("\t"+str(ind))
        file_sync = pd.read_csv(str(pop)+"_clone.vcf" ,sep="\t")
        file_sync=file_sync.sort_values(by=["#CHROM","POS"])
        for i in range(0,len(file_sync.index)):
            seq=file_sync.iloc[i]
            x=str(seq[0])
            for j in range(1,9):x=str(x)+"\t"+str(seq[j])
            for j in range(9+4*pop,9+4*pop+4):x=str(x)+"\t"+str(seq[j])
            fichier.write("\n"+str(x))
        fichier.close()
    # 3 kinship
    if kinship=="T":
        Bedgraph=open("Genotype_parent_"+str(pop)+".csv","r")
        Bedgraph=Bedgraph.read()
        Bedgraph = Bedgraph.split("\n")
        bedgraph=Bedgraph[0]
        bedgraph=bedgraph.split('\t')
        # we compare each individual
        divergence=[]
        i=4
        while i< len(bedgraph)-1:
          div=0
          tot=0
          j=i+1
          while j< len(bedgraph):
            div2=[div,tot]
            divergence.append(div2)
            j=j+1
          i=i+1
        divergence2=divergence[:]
        for line in Bedgraph[1:]: 
          divergence=divergence2[:]
          divergence2=[]
          line_plitted=line
          line_splitted = line_plitted.split("\t")
          # initial value
          ref=0
          i=4
          while i< len(bedgraph)-1:    
            j=i+1
            while j< len(bedgraph):
              x=divergence[ref]
              div=int(x[0])
              tot=int(x[1])       
              if line_splitted[i]!="NA" and line_splitted[j]!="NA":
                tot=tot+1
                if line_splitted[j]!=line_splitted[i]:div=div+1
              div2=[div,tot]
              divergence2.append(div2)
              j=j+1
              ref=ref+1
            i=i+1
        divergence=divergence2[:]
        divergence2=[]
        for div in divergence:
          div2=round(div[0]*1.0/div[1],4)
          divergence2.append(div2)
        divergence=divergence2[:]
        # we report the values
        fileout=open("Kinship_parent_"+str(pop)+".csv","a")
        fileout.write(str(len(bedgraph)-4))
        fileout.close()
        ref=0
        i=4
        while i< len(bedgraph)-1:
          written_line="\n"+str(bedgraph[i])+"\t0.0"
          j=i+1
          while j< len(bedgraph):
            x=divergence[ref]
            written_line=written_line+"\t"+str(x)
            j=j+1
            ref=ref+1
          fileout=open("Kinship_parent_"+str(pop)+".csv","a")
          fileout.write(written_line)
          fileout.close()
          i=i+1
        written_line="\n"+str(bedgraph[i])+"\t0.0"
        fileout=open("Kinship_parent_"+str(pop)+".csv","a")
        fileout.write(written_line)
        fileout.close()
        # we report the values in kinship format
        Bedgraph=open("Kinship_parent_"+str(pop)+".csv","r")
        Bedgraph=Bedgraph.read()
        Bedgraph = Bedgraph.split("\n")
        line= Bedgraph[1]
        line=line.split('\t')
        written_line=str(Bedgraph[0])+'\t'+str(line[0])
        while g<len(Bedgraph):
          line= Bedgraph[g]
          line=line.split('\t')
          written_line=written_line+'\t'+str(line[0])
          g=g+1
        fileout=open("Kinship_final_parent_"+str(pop)+".csv","a")
        fileout.write(written_line)
        fileout.close()
        fileout=open("Kinship_final_parent_"+str(pop)+".csv","a")
        line= Bedgraph[1]
        fileout.write("\n"+str(line))
        line=line.split('\t')
        divergence_ind1=line[2:]
        fileout.close()
        g=2
        while g<len(Bedgraph):
          line= Bedgraph[g]
          line=line.split('\t')
          written_line="\n"+str(line[0])
          i=1
          while i<g:
            exec("written_line=written_line+'\t'+str(divergence_ind"+str(i)+"[0])")
            exec("divergence_ind"+str(i)+"=divergence_ind"+str(i)+"[1:]")
            i=i+1
          written_line=written_line+"\t0.0"
          i=2
          while i<len(line):
            written_line=written_line+'\t'+str(line[i])
            exec("divergence_ind"+str(g)+"=line[2:]")
            i=i+1
          fileout=open("Kinship_final_parent_"+str(pop)+".csv","a")
          fileout.write(written_line)
          fileout.close()
          g=g+1
        cmd="rm Kinship_parent_"+str(pop)+".csv"
        os.system(cmd)

        Bedgraph=open("Genotype_clone_"+str(pop)+".csv","r")
        Bedgraph=Bedgraph.read()
        Bedgraph = Bedgraph.split("\n")
        bedgraph=Bedgraph[0]
        bedgraph=bedgraph.split('\t')
        # we compare each individual
        divergence=[]
        i=4
        while i< len(bedgraph)-1:
          div=0
          tot=0
          j=i+1
          while j< len(bedgraph):
            div2=[div,tot]
            divergence.append(div2)
            j=j+1
          i=i+1
        divergence2=divergence[:]
        for line in Bedgraph[1:]: 
          divergence=divergence2[:]
          divergence2=[]
          line_plitted=line
          line_splitted = line_plitted.split("\t")
          # initial value
          ref=0
          i=4
          while i< len(bedgraph)-1:    
            j=i+1
            while j< len(bedgraph):
              x=divergence[ref]
              div=int(x[0])
              tot=int(x[1])       
              if line_splitted[i]!="NA" and line_splitted[j]!="NA":
                tot=tot+1
                if line_splitted[j]!=line_splitted[i]:div=div+1
              div2=[div,tot]
              divergence2.append(div2)
              j=j+1
              ref=ref+1
            i=i+1
        divergence=divergence2[:]
        divergence2=[]
        for div in divergence:
          div2=round(div[0]*1.0/div[1],4)
          divergence2.append(div2)
        divergence=divergence2[:]
        # we report the values
        fileout=open("Kinship_clone_"+str(pop)+".csv","a")
        fileout.write(str(len(bedgraph)-4))
        fileout.close()
        ref=0
        i=4
        while i< len(bedgraph)-1:
          written_line="\n"+str(bedgraph[i])+"\t0.0"
          j=i+1
          while j< len(bedgraph):
            x=divergence[ref]
            written_line=written_line+"\t"+str(x)
            j=j+1
            ref=ref+1
          fileout=open("Kinship_clone_"+str(pop)+".csv","a")
          fileout.write(written_line)
          fileout.close()
          i=i+1
        written_line="\n"+str(bedgraph[i])+"\t0.0"
        fileout=open("Kinship_clone_"+str(pop)+".csv","a")
        fileout.write(written_line)
        fileout.close()
        # we report the values in kinship format
        Bedgraph=open("Kinship_clone_"+str(pop)+".csv","r")
        Bedgraph=Bedgraph.read()
        Bedgraph = Bedgraph.split("\n")
        line= Bedgraph[1]
        line=line.split('\t')
        written_line=str(Bedgraph[0])+'\t'+str(line[0])
        while g<len(Bedgraph):
          line= Bedgraph[g]
          line=line.split('\t')
          written_line=written_line+'\t'+str(line[0])
          g=g+1
        fileout=open("Kinship_final_clone_"+str(pop)+".csv","a")
        fileout.write(written_line)
        fileout.close()
        fileout=open("Kinship_final_clone_"+str(pop)+".csv","a")
        line= Bedgraph[1]
        fileout.write("\n"+str(line))
        line=line.split('\t')
        divergence_ind1=line[2:]
        fileout.close()
        g=2
        while g<len(Bedgraph):
          line= Bedgraph[g]
          line=line.split('\t')
          written_line="\n"+str(line[0])
          i=1
          while i<g:
            exec("written_line=written_line+'\t'+str(divergence_ind"+str(i)+"[0])")
            exec("divergence_ind"+str(i)+"=divergence_ind"+str(i)+"[1:]")
            i=i+1
          written_line=written_line+"\t0.0"
          i=2
          while i<len(line):
            written_line=written_line+'\t'+str(line[i])
            exec("divergence_ind"+str(g)+"=line[2:]")
            i=i+1
          fileout=open("Kinship_final_clone_"+str(pop)+".csv","a")
          fileout.write(written_line)
          fileout.close()
          g=g+1
        cmd="rm Kinship_clone_"+str(pop)+".csv"
        os.system(cmd)
    # 4 vcftools
    if population_analysis=="T":
        x2="dl_mean_"+str(pop)+"_parent.txt" 
        fichier3=open(str(x2), "a")
        fichier3.write("size\tR2")
        fichier3.close()
        x2="dl_mean_"+str(pop)+"_clone.txt" 
        fichier4=open(str(x2), "a")
        fichier4.write("size\tR2")
        fichier4.close()
        cmd="vcftools --vcf "+str(pop)+"_parent_sorted.vcf --geno-r2 --ld-window-bp "+str(dl_windows_size)+" --out "+str(pop)+"_dl_parent_sorted"
        os.system(cmd)
        cmd="vcftools --vcf "+str(pop)+"_clone_sorted.vcf --geno-r2 --ld-window-bp "+str(dl_windows_size)+" --out "+str(pop)+"_dl_clone_sorted"
        os.system(cmd)
        for i in range(0,(dl_windows_size)+1):
            file_sync = pd.read_csv(str(pop)+"_dl_parent_sorted.geno.ld",sep="\t")
            file_sync=file_sync[  (file_sync.POS2-file_sync.POS1==i )]
            file_sync=file_sync.dropna()
            maf=list(file_sync[file_sync.columns[4]])
            if len(maf)!=0:
              maf=sum(maf)/len(maf)
              fichier3=open(str(x2), "a")
              fichier3.write("\n"+str(i)+"\t"+str(maf))
              fichier3.close()
        for i in range(0,(dl_windows_size)+1):
            file_sync = pd.read_csv(str(pop)+"_dl_clone_sorted.geno.ld",sep="\t")
            file_sync=file_sync[  (file_sync.POS2-file_sync.POS1==i )]
            file_sync=file_sync.dropna()
            maf=list(file_sync[file_sync.columns[4]])
            if len(maf)!=0:
              maf=sum(maf)/len(maf)
              fichier4=open(str(x2), "a")
              fichier4.write("\n"+str(i)+"\t"+str(maf))
              fichier4.close()
        
        x2="pi_summary_"+str(pop)+"_parent.csv" 
        fichier=open(str(x2), "a")
        fichier.write("Pop;Chrom;pos;pi")
        x2="pi_summary_"+str(pop)+"_clone.csv" 
        fichier2=open(str(x2), "a")
        fichier2.write("Pop;Chrom;pos;pi")
        i=1
        while i <= nb_pop:
            cmd="vcftools --vcf "+str(pop)+"_parent_sorted.vcf --site-pi --keep pop"+str(i)+".txt --out "+str(pop)+"_pi_parent_sorted"+str(i)
            os.system(cmd)
            cmd="vcftools --vcf "+str(pop)+"_clone_sorted.vcf --site-pi --keep pop"+str(i)+"_clone.txt --out "+str(pop)+"_pi_clone_sorted"+str(i)
            os.system(cmd)
            cmd="vcftools --vcf "+str(pop)+"_parent_sorted.vcf --TajimaD --keep pop"+str(i)+".txt --out "+str(pop)+"_taj_parent_sorted"+str(i)
            os.system(cmd)
            cmd="vcftools --vcf "+str(pop)+"_clone_sorted.vcf --TajimaD --keep pop"+str(i)+"_clone.txt --out "+str(pop)+"_pitaj_clone_sorted"+str(i)
            os.system(cmd)
            name=str(pop)+"_pi_parent_sorted"+str(i)+".pi"
            bedgraph2 = open(str(name), "r")
            bedgraph2=bedgraph2.read()
            bedgraph2 = bedgraph2.split("\n")
            exec("pop"+str(i)+"_dict_pi_parent={}")
            for line in bedgraph2[1:]:
                if line!="":
                    line_splitted = line.split("\t")
                    pos=str(line_splitted[0])+"_"+str(line_splitted[1])
                    var=str(line_splitted[2])
                    exec("pop"+str(i)+"_dict_pi_parent[str(pos)]=var")
            name=str(pop)+"_pi_clone_sorted"+str(i)+".pi"
            bedgraph2 = open(str(name), "r")
            bedgraph2=bedgraph2.read()
            bedgraph2 = bedgraph2.split("\n")
            exec("pop"+str(i)+"_dict_pi_clone={}")
            for line in bedgraph2[1:]:
                if line!="":
                    line_splitted = line.split("\t")
                    pos=str(line_splitted[0])+"_"+str(line_splitted[1])
                    var=str(line_splitted[2])
                    exec("pop"+str(i)+"_dict_pi_clone[str(pos)]=var")
            i=i+1
        i=1
        while i <= nb_pop:
            for pos in pos_parent:
                name=str(pos).split("_")
                chrom=name[0]
                POS=name[1]
                exec("pi=pop"+str(i)+"_dict_pi_parent[str(pos)]")
                fichier.write("\npop"+str(i)+";"+str(chrom)+";"+str(POS)+";"+str(pi))
            for pos in pos_parent:
                name=str(pos).split("_")
                chrom=name[0]
                POS=name[1]
                exec("pi=pop"+str(i)+"_dict_pi_clone[str(pos)]")
                fichier2.write("\npop"+str(i)+";"+str(chrom)+";"+str(POS)+";"+str(pi))
            i=i+1
        fichier.close()
        fichier2.close()
        
        
        
        x2="fst_summary_"+str(pop)+"_parent.csv" 
        fichier5=open(str(x2), "a")
        fichier5.write("Pop;Pop2;Chrom;pos;FST")
        x2="fst_summary_"+str(pop)+"_clone.csv" 
        fichier6=open(str(x2), "a")
        fichier6.write("Pop;Pop2;Chrom;pos;FST")
        i=1
        while i <= (nb_pop-1):
            j=i+1
            while j <= (nb_pop):
                cmd="vcftools --vcf "+str(pop)+"_parent_sorted.vcf --weir-fst-pop pop"+str(i)+".txt --weir-fst-pop pop"+str(j)+".txt --out "+str(pop)+"_parent_FST_pop"+str(i)+"_pop"+str(j)+"_geno"
                os.system(cmd)
                cmd="vcftools --vcf "+str(pop)+"_clone_sorted.vcf --weir-fst-pop pop"+str(i)+".txt --weir-fst-pop pop"+str(j)+".txt --out "+str(pop)+"_clone_FST_pop"+str(i)+"_pop"+str(j)+"_geno"
                os.system(cmd)
                name=str(pop)+"_parent_FST_pop"+str(i)+"_pop"+str(j)+"_geno.fst"
                bedgraph2 = open(str(name), "r")
                bedgraph2=bedgraph2.read()
                bedgraph2 = bedgraph2.split("\n")
                exec("pop"+str(i)+"_pop"+str(j)+"_dict_fst_parent={}")
                for line in bedgraph2[1:]:
                    if line!="":
                        line_splitted = line.split("\t")
                        pos=str(line_splitted[0])+"_"+str(line_splitted[1])
                        var=str(line_splitted[2])
                        exec("pop"+str(i)+"_pop"+str(j)+"_dict_fst_parent[str(pos)]=var")
                name=str(pop)+"_clone_FST_pop"+str(i)+"_pop"+str(j)+"_geno.fst"
                bedgraph2 = open(str(name), "r")
                bedgraph2=bedgraph2.read()
                bedgraph2 = bedgraph2.split("\n")
                exec("pop"+str(i)+"_pop"+str(j)+"_dict_fst_clone={}")
                for line in bedgraph2[1:]:
                    if line!="":
                        line_splitted = line.split("\t")
                        pos=str(line_splitted[0])+"_"+str(line_splitted[1])
                        var=str(line_splitted[2])
                        exec("pop"+str(i)+"_pop"+str(j)+"_dict_fst_clone[str(pos)]=var")
                j=j+1
            i=i+1
            
        i=1
        while i <= (nb_pop-1):
            j=i+1
            while j <= (nb_pop):
                for pos in pos_parent:
                    name=str(pos).split("_")
                    chrom=name[0]
                    POS=name[1]
                    exec("pi=pop"+str(i)+"_pop"+str(j)+"_dict_fst_parent[str(pos)]")
                    fichier5.write("\npop"+str(i)+";pop"+str(j)+";"+str(chrom)+";"+str(POS)+";"+str(pi))
                    exec("pi=pop"+str(i)+"_pop"+str(j)+"_dict_fst_clone[str(pos)]")
                    fichier6.write("\npop"+str(i)+";pop"+str(j)+";"+str(chrom)+";"+str(POS)+";"+str(pi))
                j=j+1
            i=i+1
        fichier5.close()
        fichier6.close()
    # conservation
    if conservation=="T":
        x2=str(pop)+"_conservation_ind.csv" 
        output=open(str(x2), "a")
        output.write("ind\ttot\tnb_conserv_met\tnb_conserv_demet\tnb_ant_met\tnb_ant_demet")
        filein1=open("Genotype_parent_"+str(pop)+".csv","r")
        filein1=filein1.read()
        filein1=filein1.split('\n')
        filein1=filein1[1:]
        filein2=open("Genotype_clone_"+str(pop)+".csv","r")
        filein2=filein2.read()
        filein2=filein2.split('\n')
        filein2=filein2[1:]
        for ind in range(4,(4+len(id_ind))):
            ID=id_ind[ind-4]
            pos_parent={}
            for line in filein1:
                x=line.split("\t")
                if x[i]!="NA":
                    pos=str(x[2])+"_"+str(x[3])
                    pos_parent[pos]=str(x[i])
            met=0
            demet=0
            tot=0
            ant_met=0
            ant_dem=0
            for line in filein2:
                x=line.split("\t")
                pos=str(x[2])+"_"+str(x[3])
                if x[i]!="NA" and pos in list(pos_parent.keys()):
                    tot=tot+1
                    x2=pos_parent[pos]
                    if str(x[i])!=x2 and x2=="A":ant_met=ant_met+1
                    elif str(x[i])!=x2 and x2=="T":ant_dem=ant_dem+1
                    elif str(x[i])==x2 and x2=="A":met=met+1
                    elif str(x[i])==x2 and x2=="T":demet=demet+1


            output.write("\n"+str(ID)+"\t"+str(tot)+"\t"+str(met)+"\t"+str(demet)+"\t"+str(ant_met)+"\t"+str(ant_dem))
            



        
                        
            
        
    

    
    
if __name__=='__main__':
    
    ArgsVal = sys.argv[1:]
    Diversity_methylation_program(ArgsVal)
