VCF="../data/vcf_files"
VCF_MERGED="../data/vcf_files_control_treatment_merged"
MAF_MERGED="../data/maf_vep_files_control_treatment_merged"
S_LIST="../data/sample_metadata"
VCF2MAF="/home/dell/mskcc-vcf2maf-958809e/vcf2maf.pl"
REMAP="/home/dell/mskcc-vcf2maf-958809e/data/hg19_to_GRCh37.chain"

from glob import glob
import os
import os.path
from os import path

class GenerateMaf(object):
    
    def __init__(self):
        """Initialization."""
        pass

    def run_cell_line_analysis(self,VCF,VCF_MERGED,MAF_MERGED):
        subdir_list=glob("%s/*/"%VCF)
        for sub_directory in subdir_list:
            dir_suffix=self.make_merged_subdirectory(VCF_MERGED,sub_directory)
            metadata_dictionary=self.make_sample_metadata_dictionary(S_LIST,dir_suffix)
            self.change_vcf_ids(metadata_dictionary,dir_suffix)
            self.merge_bcftools_vcf_condition_control(metadata_dictionary,dir_suffix)
            self.vcf_to_maf(metadata_dictionary,dir_suffix)

    def vcf_to_maf(self,metadata_dictionary,dir_suffix):
        os.system("mkdir -p %s/%s_MAF"%(MAF_MERGED,dir_suffix))
        for filename,metaid in metadata_dictionary.items():
            if not metaid.startswith("control"):
                tumor_id=metaid;control_id="control_"+"_".join(metaid.split("_")[1:])
                vcf2maf="perl %s --input-vcf %s/%s_VCF/%s_merged_control_%s.vcf --output-maf %s/%s_MAF/%s_merged_control_%s.vcf.maf --remap-chain %s --tumor-id %s --vcf-tumor-id %s --normal-id %s --vcf-normal-id %s"%(VCF2MAF,VCF_MERGED,dir_suffix,dir_suffix,metaid,MAF_MERGED,dir_suffix,dir_suffix,metaid,REMAP,tumor_id,tumor_id,control_id,control_id)
                os.system(vcf2maf)

    def merge_bcftools_vcf_condition_control(self,metadata_dictionary,dir_suffix):
        reverse_dictionary=dict((v,k) for k,v in metadata_dictionary.items())
        os.system("mkdir -p %s/%s_VCF"%(VCF_MERGED,dir_suffix))

        for metaid,filename in reverse_dictionary.items():
            if metaid.startswith("control"):
                rep_id="_".join(metaid.split("_")[1:])
                control_metaid="control_"+rep_id; control_filename=reverse_dictionary[control_metaid]
                id_changed_control_file=control_filename[:-4]+"_id_changed.vcf"
                
                os.system("bgzip %s/%s_VCF/%s"%(VCF,dir_suffix,id_changed_control_file))
                os.system("bcftools index %s/%s_VCF/%s.gz"%(VCF,dir_suffix,id_changed_control_file))

        for metaid,filename in reverse_dictionary.items():
            if not metaid.startswith("control"):
                id_changed_condition_file=filename[:-4]+"_id_changed.vcf"

                rep_id="_".join(metaid.split("_")[1:])
                control_metaid="control_"+rep_id; control_filename=reverse_dictionary[control_metaid]

                id_changed_control_file=control_filename[:-4]+"_id_changed.vcf"


                os.system("bgzip %s/%s_VCF/%s"%(VCF,dir_suffix,id_changed_condition_file))
                os.system("bcftools index %s/%s_VCF/%s.gz"%(VCF,dir_suffix,id_changed_condition_file))

                os.system("bcftools merge %s/%s_VCF/%s.gz %s/%s_VCF/%s.gz -Ov -o %s/%s_VCF/%s_merged_control_%s.vcf"%(VCF,dir_suffix,id_changed_control_file,VCF,dir_suffix,id_changed_condition_file,VCF_MERGED,dir_suffix,dir_suffix,metaid))

    def change_vcf_ids(self,metadata_dictionary,dir_suffix):
        for filename,metaid in metadata_dictionary.items():
            file_content=[]
            file_content_modified=[]

            if not path.exists("%s/%s_VCF/%s"%(VCF,dir_suffix,filename)):
                os.system("gunzip -k %s/%s_VCF/%s.gz"%(VCF,dir_suffix,filename))

            with open("%s/%s_VCF/%s"%(VCF,dir_suffix,filename)) as fh:
                file_content=fh.readlines()
                for line in file_content:
                    if line.startswith("#CHROM\t"):
                        line_list=line.strip().split("\t")
                        line_list[-1]=metaid
                        line_str="\t".join(line_list)
                        file_content_modified.append(line_str)

                    else:file_content_modified.append(line.strip())

            with open("%s/%s_VCF/%s_id_changed.vcf"%(VCF,dir_suffix,filename[:-4]),"w") as fh:
                for line in file_content_modified:
                    fh.write("%s\n"%line)

    def make_sample_metadata_dictionary(self,S_LIST,dir_suffix):
        metadata_dictionary={}
        with open("%s/%s_sample_list.tsv"%(S_LIST,dir_suffix.lower())) as fh:
            next(fh)
            for line in fh:
                line_list=line.rstrip().split("\t")
                filename=line_list[0]
                metaid=line_list[1:]
                metadata_dictionary[filename]="_".join(metaid)
        return metadata_dictionary

    def make_merged_subdirectory(self,VCF_MERGED,sub_directory):
        dir_suffix=sub_directory.split("/")[-2].split("_")[0]
        mkdir="mkdir -p %s/%s_VCF"%(VCF_MERGED,dir_suffix)
        os.system(mkdir)
        return dir_suffix

if __name__=="__main__":
    run_pipeline=GenerateMaf()
    run_pipeline.run_cell_line_analysis(VCF,VCF_MERGED,MAF_MERGED)
