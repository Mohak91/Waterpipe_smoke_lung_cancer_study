#input files
WD="/media/dell/MyPassport/Somatic_Cancer_analysis_wo_merging_technical_replicate_bams_17_may_2021/scripts/TCGA_tumor_mutational_burden_analysis"
BOTH_IDS="data_bcr_clinical_data_sample.txt"
MAF="data_mutations_mskcc.txt"
SHORTER_ID="data_bcr_clinical_data_patient.txt" #also contains meta data like age, gender, tumor stage etc.

#files generated from the script
SMOKING_MAF="data_mutations_mskcc_smoking.maf"
NON_SMOKING_MAF="data_mutations_mskcc_non_smoking.maf"
R_SCRIPT="tumor_mutational_burden_tcga_smoking_non_smoking.R"
TMB_CSV="tumor_mutational_burden_smoker_non_smokers.csv"
TMB_BOXPLOT_PDF="tumor_mutational_burden_boxplot.pdf"
INPUT_FOR_RFFI="tumor_mutational_burden_random_forest_feature_importance_analysis_input_table.tsv"
R_SCRIPT_RF="tmb_random_forest_feature_importance_analysis.R"
TMB_RF_FEATURE_IMPORTANCE_PDF="percentage_incMSE_tmb_factors.pdf"

__prog_name__ = 'tumor_mutational_burden_analysis_tcga_between_cohorts.py'
__prog_desc__ = 'Compare Tumor Mutational Burden on a TCGA dataset between smokers and non-smokers and run random forest regression based feature importance for patient metadata.'

__author__ = 'Mohak Sharda'
__copyright__ = 'Copyright 2022'
__credits__ = ['Mohak Sharda']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Mohak Sharda'
__email__ = 'mohaks@ncbs.res.in'
__status__ = 'Development'

import os

class TumorMutationalBurdenAnalysis(object):

    patient_sample_dictionary_class={}
    patient_id_other_confounding_factors_dictionary_class={}

    def run_tmb_analysis(self, maf_file, sample_file, patient_file):
        patient_sample_dictionary=self.map_patient_sample_identifiers(sample_file)
        patient_id_smoking_status_dictionary,patient_id_other_confounding_factors_dictionary=self.map_patient_id_smoking_status(patient_file)
        TumorMutationalBurdenAnalysis.patient_sample_dictionary_class=patient_sample_dictionary
        TumorMutationalBurdenAnalysis.patient_id_other_confounding_factors_dictionary_class=patient_id_other_confounding_factors_dictionary
        header,sample_id_maf_dictionary=self.map_sample_id_maf_rows(maf_file)
        sample_id_smoking_status_dictionary=self.merge_sample_id_smoking_status(patient_sample_dictionary,patient_id_smoking_status_dictionary)
        smoking_maf,non_smoking_maf=self.segregate_dictionaries_smoking_status(sample_id_smoking_status_dictionary,sample_id_maf_dictionary)
        self.write_file(header,smoking_maf,SMOKING_MAF)
        self.write_file(header,non_smoking_maf,NON_SMOKING_MAF)
        self.prepare_r_script(SMOKING_MAF,NON_SMOKING_MAF,R_SCRIPT,TMB_CSV)
        self.run_r_script(R_SCRIPT)

    def map_patient_sample_identifiers(self, sample_file):
        patient_sample_dictionary={}
        with open(sample_file) as fh:
            for line in fh:
                if not line.startswith("#") and not line.startswith("PATIENT_ID"):
                    line_list=line.rstrip().split("\t")
                    patient_id=line_list[0]
                    sample_id=line_list[1]
                    patient_sample_dictionary[patient_id]=sample_id

        return patient_sample_dictionary

    def map_patient_id_smoking_status(self, patient_file):
        patient_id_smoking_status_dictionary={}
        patient_id_other_confounding_factors_dictionary={}
        with open(patient_file) as fh:
            for line in fh:
                if not line.startswith("#") and not line.startswith("OTHER_PATIENT_ID"):
                    line_list=line.rstrip().split("\t")
                    patient_id=line_list[1]
                    smoking_status=line_list[40]
                    gender=line_list[6]
                    age=line_list[53]
                    ajcc_staging=line_list[19]; ajcc_nodes=line_list[20];metastasis=line_list[21]; ajcc_pathology=line_list[22]
                    patient_id_smoking_status_dictionary[patient_id]=smoking_status
                    patient_id_other_confounding_factors_dictionary[patient_id]=[smoking_status,gender,age,ajcc_staging,ajcc_nodes,metastasis,ajcc_pathology]

        return patient_id_smoking_status_dictionary,patient_id_other_confounding_factors_dictionary
    
    def map_sample_id_maf_rows(self,maf_file):
        sample_id_maf_dictionary={};header=[]
        with open(maf_file) as fh:
            header=fh.readline().rstrip().split("\t")
            for line in fh:
                line_list=line.rstrip().split("\t")
                sample_id=line_list[16]
                if not sample_id in sample_id_maf_dictionary:
                    sample_id_maf_dictionary[sample_id]=[line_list]
                elif sample_id in sample_id_maf_dictionary:
                    sample_id_maf_dictionary[sample_id].append(line_list)
        
        return header,sample_id_maf_dictionary

    def merge_sample_id_smoking_status(self,patient_sample_dictionary,patient_id_smoking_status_dictionary):
        sample_id_smoking_status_dictionary={}
        for patient_id,sample_id in patient_sample_dictionary.items():
            if patient_id in patient_id_smoking_status_dictionary.keys():
                smoking_status=patient_id_smoking_status_dictionary[patient_id]
                sample_id_smoking_status_dictionary[sample_id]=smoking_status

        return sample_id_smoking_status_dictionary

    def segregate_dictionaries_smoking_status(self,sample_id_smoking_status_dictionary,sample_id_maf_dictionary):
        smoking_maf=[]
        non_smoking_maf=[]
        for sample_id,smoking_status in sample_id_smoking_status_dictionary.items():
            if smoking_status=="1" and sample_id in sample_id_maf_dictionary:
                for maf_rows in sample_id_maf_dictionary[sample_id]:
                    non_smoking_maf.append(maf_rows)
            elif smoking_status=="2" and sample_id in sample_id_maf_dictionary:
                for maf_rows in sample_id_maf_dictionary[sample_id]:
                    smoking_maf.append(maf_rows)

        return smoking_maf,non_smoking_maf

    def write_file(self,header,to_write_list,write_file):
        with open(write_file,"w") as fh:
            header_string="\t".join(str(x) for x in header)
            fh.write("%s\n"%header_string)
            for maf_row_list in to_write_list:
                maf_row_string="\t".join(str(x) for x in maf_row_list)
                fh.write("%s\n"%maf_row_string)

    def prepare_r_script(self,SMOKING_MAF,NON_SMOKING_MAF,R_SCRIPT,TMB_CSV):
        with open(R_SCRIPT,"w") as fh:
            fh.write("library(maftools)\n")
            fh.write("setwd(\"%s\")\n"%WD)
            fh.write("smoke_maf<-read.maf(c(\"%s\"))\n"%SMOKING_MAF)
            fh.write("non_smoke_maf<-read.maf(c(\"%s\"))\n"%NON_SMOKING_MAF)
            fh.write("mutload = tcgaCompare(maf = c(non_smoke_maf,smoke_maf), cohortName = c(\"non_smokers\",\"smokers\"), logscale = TRUE, capture_size = 50)\n")
            fh.write("write.csv(x=mutload$mutation_burden_perSample,file = \"%s\")\n"%TMB_CSV)
            fh.write("mutload_read<-read.table(\"%s\",sep=\",\",header=T)\n"%TMB_CSV)
            fh.write("library(dplyr)\n")
            fh.write("nonsmokers_smokers_cohorts<-mutload_read %>% filter(cohort==\"non_smokers\" | cohort==\"smokers\")\n")
            fh.write("plot_non_smokers_smokers_cohorts <- nonsmokers_smokers_cohorts %>% ggplot(aes(x=cohort,y=total_perMB,fill=factor(cohort))) + geom_boxplot(width=0.5,lwd=1) + geom_jitter(width=0.15)+ labs(subtitle=\"TMB (per MB)\")\n")
            fh.write("pdf(\"%s\")\n"%TMB_BOXPLOT_PDF)
            fh.write("plot_non_smokers_smokers_cohorts + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour=\"black\"))\n")
            fh.write("dev.off()\n")
            fh.write("group_by(nonsmokers_smokers_cohorts, cohort) %>% summarise(count = n(), median = median(total_perMB, na.rm = TRUE), IQR = IQR(total_perMB, na.rm = TRUE))\n")
            fh.write("wilcox.test(total_perMB ~ cohort, data = nonsmokers_smokers_cohorts,exact = FALSE)\n")

    def run_r_script(self,R_SCRIPT):
        os.system("Rscript %s"%R_SCRIPT)

class RandomForestFeatureImportanceAnalysis(TumorMutationalBurdenAnalysis):

    def run_rffi_analysis(self,tmb_csv_file):
        sample_id_tmb_dictionary=self.extract_sample_id_tmb_from_csv_file(tmb_csv_file)
        sample_id_other_confounding_factors_dictionary=self.merge_sample_id_meta_data(TumorMutationalBurdenAnalysis.patient_sample_dictionary_class,TumorMutationalBurdenAnalysis.patient_id_other_confounding_factors_dictionary_class)
        self.write_to_file(sample_id_tmb_dictionary,sample_id_other_confounding_factors_dictionary,INPUT_FOR_RFFI)
        self.prepare_random_forest_feature_importance_R_script(INPUT_FOR_RFFI,R_SCRIPT_RF)
        self.run_r_script(R_SCRIPT_RF)

    def merge_sample_id_meta_data(self,patient_sample_dictionary,patient_id_other_confounding_factors_dictionary):
        sample_id_other_confounding_factors_dictionary={}
        for patient_id,sample_id in patient_sample_dictionary.items():
            if patient_id in patient_id_other_confounding_factors_dictionary.keys():
                all_confounding_factors=patient_id_other_confounding_factors_dictionary[patient_id] #includes smoking status
                sample_id_other_confounding_factors_dictionary[sample_id]=all_confounding_factors #value is a list

        return sample_id_other_confounding_factors_dictionary

    def extract_sample_id_tmb_from_csv_file(self,tmb_csv_file):
        sample_id_tmb_dictionary={}
        with open(tmb_csv_file) as fh:
            next(fh)
            for line in fh:
                line_list=line.rstrip().split(",")
                sample_id=line_list[1][1:-1]
                tmb=line_list[4]
                sample_id_tmb_dictionary[sample_id]=tmb

        return sample_id_tmb_dictionary

    def write_to_file(self,sample_id_tmb_dictionary,sample_id_other_confounding_factors_dictionary,to_write_file_for_rffi):
        with open(to_write_file_for_rffi,"w") as fh:
            fh.write("TMB\tsmoking_status\tgender\tage\tajcc_staging\tajcc_nodes\tmetastasis\tajcc_pathology\n")
            for sample_id,all_confounding_factors in sample_id_other_confounding_factors_dictionary.items():
                if sample_id in sample_id_tmb_dictionary:
                    tmb=sample_id_tmb_dictionary[sample_id]
                    other_factors_string="\t".join(str(x) for x in all_confounding_factors)
                    if not "[" in other_factors_string and not "]" in other_factors_string: #filters [non-available] meta data rows
                        fh.write("%s\t%s\n"%(tmb,other_factors_string))

    def prepare_random_forest_feature_importance_R_script(self,INPUT_FOR_RFFI,R_SCRIPT_RF):
        with open(R_SCRIPT_RF,"w") as fh:
            fh.write("library(dplyr)\n")
            fh.write("library(ggplot2)\n")
            fh.write("library(randomForest)\n")
            #Create the dataframe
            fh.write("tmb_table <- as.data.frame(read.table(\"%s\",sep=\"\\t\",header=T))\n"%INPUT_FOR_RFFI)
            fh.write("tmb_table[tmb_table==\"T1\"]<-\"0\"\n")
            fh.write("tmb_table[tmb_table==\"T1a\"]<-\"1\"\n")
            fh.write("tmb_table[tmb_table==\"T1b\"]<-\"2\"\n")
            fh.write("tmb_table[tmb_table==\"T2\"]<-\"3\"\n")
            fh.write("tmb_table[tmb_table==\"T2a\"]<-\"4\"\n")
            fh.write("tmb_table[tmb_table==\"T2b\"]<-\"5\"\n")
            fh.write("tmb_table[tmb_table==\"T3\"]<-\"6\"\n")
            fh.write("tmb_table[tmb_table==\"T4\"]<-\"7\"\n")
            fh.write("tmb_table[tmb_table==\"NX\"]<-\"0\"\n")
            fh.write("tmb_table[tmb_table==\"N0\"]<-\"1\"\n")
            fh.write("tmb_table[tmb_table==\"N1\"]<-\"2\"\n")
            fh.write("tmb_table[tmb_table==\"N2\"]<-\"3\"\n")
            fh.write("tmb_table[tmb_table==\"N3\"]<-\"4\"\n")
            fh.write("tmb_table[tmb_table==\"MX\"]<-\"0\"\n")
            fh.write("tmb_table[tmb_table==\"M0\"]<-\"1\"\n")
            fh.write("tmb_table[tmb_table==\"M1\"]<-\"2\"\n")
            fh.write("tmb_table[tmb_table==\"Stage I\"]<-\"0\"\n")
            fh.write("tmb_table[tmb_table==\"Stage IA\"]<-\"1\"\n")
            fh.write("tmb_table[tmb_table==\"Stage IB\"]<-\"2\"\n")
            fh.write("tmb_table[tmb_table==\"Stage IIA\"]<-\"3\"\n")
            fh.write("tmb_table[tmb_table==\"Stage IIB\"]<-\"4\"\n")
            fh.write("tmb_table[tmb_table==\"Stage IIIA\"]<-\"5\"\n")
            fh.write("tmb_table[tmb_table==\"Stage IIIB\"]<-\"6\"\n")
            fh.write("tmb_table[tmb_table==\"Stage IV\"]<-\"7\"\n")
            fh.write("tmb_table <- transform(tmb_table,TMB=as.numeric(TMB),smoking_status=as.factor(smoking_status),gender=as.factor(gender),age=as.numeric(age),ajcc_staging=as.ordered(ajcc_staging),ajcc_nodes=as.ordered(ajcc_nodes),metastasis=as.ordered(metastasis),ajcc_pathology=as.ordered(ajcc_pathology))\n")
            fh.write("model <- randomForest(TMB ~ smoking_status+gender+age+ajcc_staging+ajcc_nodes+metastasis+ajcc_pathology,data = tmb_table, importance=TRUE)\n")
            fh.write("colnames(model$importance)<-c(\"IncMSE\",\"IncNodePurity\")\n")
            fh.write("dataframe_model<-as.data.frame(model$importance)\n")
            fh.write("dataframe_model$Factors<-row.names(dataframe_model)\n")
            fh.write("pdf(\"%s\")\n"%TMB_RF_FEATURE_IMPORTANCE_PDF)
            fh.write("ggplot(dataframe_model,aes(x=Factors,y=IncMSE))+geom_point()+geom_segment(aes(x=Factors,xend=Factors,y=0,yend=IncMSE))\n")
            fh.write("dev.off()\n")

    def run_r_script(self,R_SCRIPT_RF):
        os.system("Rscript %s"%R_SCRIPT_RF)

if __name__=="__main__":
    
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')
    
    tcga_dataset=TumorMutationalBurdenAnalysis()
    tcga_dataset.run_tmb_analysis(MAF,BOTH_IDS,SHORTER_ID)
    tmb_dataset=RandomForestFeatureImportanceAnalysis()
    tmb_dataset.run_rffi_analysis(TMB_CSV)
