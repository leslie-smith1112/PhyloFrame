#include <iostream>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>
#include <hts.h>
#include <vcf.h>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

//Spacing and open/close functions from htslib example: http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
//Variable naimg from samtools

int main(int argc, char **argv) {

    CLI::App app("Frequencies reader");
    std::string vcf_file_name;
    std::string out_file;
    app.add_option("-i,--input", vcf_file_name, "Input VCF file")->required();
    app.add_option("-o,--output", out_file, "Output file")->required();
    CLI11_PARSE(app, argc, argv);

    htsFile* inf = bcf_open(vcf_file_name.c_str(), "r");
    if (inf == NULL) {
        spdlog::error("Error in open of vcf file: {}", vcf_file_name);
        std::exit(EXIT_FAILURE);
    }

    spdlog::info("Parsing vcf: {}", vcf_file_name);

    bcf_hdr_t* hdr = bcf_hdr_read(inf);
    bcf1_t* rec = bcf_init();
    if (rec == NULL) {
        spdlog::error("Error parsing vcf file: {}", vcf_file_name);
        bcf_close(inf);
        bcf_hdr_destroy(hdr);
        std::exit(EXIT_FAILURE);
    }


    // gnomad4.1 exomes does not have amish. The genomes do.
    // Removed from below: "AMI_AC" << "\t" << "AMI_AN" << "\t" << "AMI_AF" << "\t" << 
    std::ofstream out_file_stream(out_file);
    if (not out_file_stream.is_open()) { spdlog::error("Can't open output file {}", out_file); }
    else { spdlog::info("Writing output to {}", out_file); }
    //out_file_stream << "CHROM" << "\t" << "POS" << "\t" << "RSID" <<"\t" << "REF" << "\t" << "ALT" << "\t" << "AC" << "\t" <<
                    //"AN" << "\t" << "AF" << "\t" << "EAS_AF" << "\t" << "AMR_AF" << "\t" << "AFR_AF" << "\t" << "EUR_AF" << "\t" << "SAS_AF" << "\n";
    out_file_stream << "CHROM" << "\t" << "POS" << "\t" << "RSID" <<"\t" << "REF" << "\t" << "ALT" << "\t" <<
    "AFR_AC" << "\t" << "AFR_AN" << "\t" << "AFR_AF" << "\t" << 
    "AMR_AC" << "\t" << "AMR_AN" << "\t" << "AMR_AF" << "\t" << 
    "ASJ_AC"  << "\t"<< "ASJ_AN" << "\t" << "ASJ_AF" << "\t" << 
    "EAS_AC" << "\t" << "EAS_AN" << "\t" << "EAS_AF" << "\t" << 
    "FIN_AC" << "\t" << "FIN_AN" << "\t" << "FIN_AF" << "\t" << 
    "MID_AC" << "\t" << "MID_AN" << "\t" << "MID_AF" << "\t" << 
    "NFE_AC" << "\t" << "NFE_AN" << "\t" << "NFE_AF" << "\t" << 
    "SAS_AC" << "\t" << "SAS_AN" << "\t" << "SAS_AF" << "\t" << 
    "CADD_PHRED" << "\t" << "REVEL_MAX" << "\t" << "SPLICEAI_DS_MAX" << 
    "\t" << "PANGOLIN_LARGEST_DS" << "\t" << "PHYLOP" << "\t" << "SIFT_MAX" << 
    "\t" << "POLYPHEN_MAX" << "\t" << "VEP_VAL" << "\n";

   // start parsing
    while (bcf_read(inf, hdr, rec) == 0) {
        //cout << "here";

        std::string line;
        bcf_unpack(rec, BCF_UN_ALL);
        int chrom = rec->rid;
        int position = rec->pos;
        //cout << "position: " << position;
        const char* chr_name = bcf_hdr_id2name(hdr, rec->rid);
//      bcf_unpack(rec,BCF_UN_STR);
        bcf_dec_t qhat = rec ->d;
        char* ID = qhat.id;
        std::string someString(ID);
        out_file_stream << chr_name << "\t"<< position << "\t"<< ID << "\t";
        //get reference and alternate allele
        std::string ref_allele = rec->d.allele[0];
        std::string alt_allele = rec->d.allele[1];
        bcf_info_t* info_pointer  = bcf_get_info(hdr, rec, "vep");
        std::string vep_val((char*)info_pointer->vptr);
        out_file_stream << ref_allele << "\t" << alt_allele << "\t";

   //     bcf_unpack(rec, BCF_UN_ALL);

        //GNOMADV4 VERSION::
        //AFRICAN
        bcf_info_t* AN_afr  = bcf_get_info(hdr, rec, "AN_afr");
        bcf_info_t* AC_afr  = bcf_get_info(hdr, rec, "AC_afr");
        bcf_info_t* AF_afr  = bcf_get_info(hdr, rec, "AF_afr");
        int AC_afr_2 = AC_afr->v1.i; // allele count will always be an integer
        int AN_afr_2 = AN_afr->v1.i; // allele nufber will always be an interger
        float AF_afr_2;
        if(AF_afr == NULL) {
            if(AN_afr_2 != 0){
                AF_afr_2 = AC_afr_2 / AN_afr_2;
            }
            else{
                AF_afr_2 = 0;
            }
        }else{
            AF_afr_2 = AF_afr->v1.f;
        }
        out_file_stream <<  AC_afr_2 << "\t";
        out_file_stream <<  AN_afr_2 << "\t";
        out_file_stream <<  AF_afr_2 << "\t";

        //LATINO
        bcf_info_t* AN_amr  = bcf_get_info(hdr, rec, "AN_amr");
        bcf_info_t* AC_amr  = bcf_get_info(hdr, rec, "AC_amr");
        bcf_info_t* AF_amr  = bcf_get_info(hdr, rec, "AF_amr");
        int AC_amr_2 = AC_amr->v1.i; // allele count will always be an integer
        int AN_amr_2 = AN_amr->v1.i; // allele number will always be an interger

        //assign AF to the AF value in the file if it exists, otherwise try to calculate it
        float AF_amr_2;  // allele frequency should always be a float (I dont understand what the integer value returns here ?? )
        if(AF_amr == NULL) {
            if(AN_amr_2 != 0){
                AF_amr_2 = AC_amr_2 / AN_amr_2;
            }
            else{
                AF_amr_2 = 0;
            }
        }else{
            AF_amr_2 = AF_amr->v1.f;
        }

        out_file_stream <<  AC_amr_2 << "\t";
        out_file_stream <<  AN_amr_2 << "\t";
        out_file_stream <<  AF_amr_2 << "\t";

        //excluded for gnomad exomes 4.1 
        /*
        //AMISH
        bcf_info_t* AN_ami  = bcf_get_info(hdr, rec, "AN_ami");
        bcf_info_t* AC_ami  = bcf_get_info(hdr, rec, "AC_ami");
        bcf_info_t* AF_ami  = bcf_get_info(hdr, rec, "AF_ami");
        int AC_ami_2 = AC_ami->v1.i; // allele count will always be an integer
        int AN_ami_2 = AN_ami->v1.i; // allele nufber will always be an interger

        float AF_ami_2;
        if(AF_ami == NULL) {
            if(AN_ami_2 != 0){
                AF_ami_2 = AC_ami_2 / AN_ami_2;
            }
            else{
                AF_ami_2 = 0;
            }
        }else{
            AF_ami_2 = AF_ami->v1.f;
        }

        out_file_stream <<  AC_ami_2 << "\t";
        out_file_stream <<  AN_ami_2 << "\t";
        out_file_stream <<  AF_ami_2 << "\t";
        */

        // ASHKENAZI JEWISH
        bcf_info_t* AN_asj  = bcf_get_info(hdr, rec, "AN_asj");
        bcf_info_t* AC_asj  = bcf_get_info(hdr, rec, "AC_asj");
        bcf_info_t* AF_asj  = bcf_get_info(hdr, rec, "AF_asj");
        int AC_asj_2 = AC_asj->v1.i; // allele count will always be an integer
        int AN_asj_2 = AN_asj->v1.i; // allele nufber will always be an interger

        float AF_asj_2;
        if(AF_asj == NULL) {
            if(AN_asj_2 != 0){
                AF_asj_2 = AC_asj_2 / AN_asj_2;
            }
            else{
                AF_asj_2 = 0;
            }
        }else{
            AF_asj_2 = AF_asj->v1.f;
        }

        out_file_stream <<  AC_asj_2 << "\t";
        out_file_stream <<  AN_asj_2 << "\t";
        out_file_stream <<  AF_asj_2 << "\t";


        //EAST ASIAN
        bcf_info_t* AN_eas  = bcf_get_info(hdr, rec, "AN_eas");
        bcf_info_t* AC_eas  = bcf_get_info(hdr, rec, "AC_eas");
        bcf_info_t* AF_eas  = bcf_get_info(hdr, rec, "AF_eas");
        int AC_eas_2 = AC_eas->v1.i; // allele count will always be an integer
        int AN_eas_2 = AN_eas->v1.i; // allele nufber will always be an interger

        float AF_eas_2;
        if(AF_eas == NULL) {
            if(AN_eas_2 != 0){
                AF_eas_2 = AC_eas_2 / AN_eas_2;
            }
            else{
                AF_eas_2 = 0;
            }
        }else{
            AF_eas_2 = AF_eas->v1.f;
        }

        out_file_stream <<  AC_eas_2 << "\t";
        out_file_stream <<  AN_eas_2 << "\t";
        out_file_stream <<  AF_eas_2 << "\t";

        //FINNISH
        bcf_info_t* AN_fin  = bcf_get_info(hdr, rec, "AN_fin");
        bcf_info_t* AC_fin  = bcf_get_info(hdr, rec, "AC_fin");
        bcf_info_t* AF_fin  = bcf_get_info(hdr, rec, "AF_fin");
        int AC_fin_2 = AC_fin->v1.i; // allele count will always be an integer
        int AN_fin_2 = AN_fin->v1.i; // allele nufber will always be an interger

        float AF_fin_2;
        if(AF_fin == NULL) {
            if(AN_fin_2 != 0){
                AF_fin_2 = AC_fin_2 / AN_fin_2;
            }
            else{
                AF_fin_2 = 0;
            }
        }else{
            AF_fin_2 = AF_fin->v1.f;
        }

        out_file_stream <<  AC_fin_2 << "\t";
        out_file_stream <<  AN_fin_2 << "\t";
        out_file_stream <<  AF_fin_2 << "\t";

        //MIDDLE EASTERN
        bcf_info_t* AN_mid  = bcf_get_info(hdr, rec, "AN_mid");
        bcf_info_t* AC_mid  = bcf_get_info(hdr, rec, "AC_mid");
        bcf_info_t* AF_mid  = bcf_get_info(hdr, rec, "AF_mid");
        int AC_mid_2 = AC_mid->v1.i; // allele count will always be an integer
        int AN_mid_2 = AN_mid->v1.i; // allele nufber will always be an interger

        float AF_mid_2;
        if(AF_mid == NULL) {
            if(AN_mid_2 != 0){
                AF_mid_2 = AC_mid_2 / AN_mid_2;
            }
            else{
                AF_mid_2 = 0;
            }
        }else{
            AF_mid_2 = AF_mid->v1.f;
        }

        out_file_stream <<  AC_mid_2 <<"\t";
        out_file_stream <<  AN_mid_2 <<"\t";
        out_file_stream <<  AF_mid_2 << "\t";

        //NONFINNISH EUROPEAN
        bcf_info_t* AN_nfe  = bcf_get_info(hdr, rec, "AN_nfe");
        bcf_info_t* AC_nfe  = bcf_get_info(hdr, rec, "AC_nfe");
        bcf_info_t* AF_nfe  = bcf_get_info(hdr, rec, "AF_nfe");
        int AC_nfe_2 = AC_nfe->v1.i; // allele count will always be an integer
        int AN_nfe_2 = AN_nfe->v1.i; // allele nufber will always be an interger

        float AF_nfe_2;
        if(AF_nfe == NULL) {
            if(AN_nfe_2 != 0){
                AF_nfe_2 = AC_nfe_2 / AN_nfe_2;
            }
            else{
                AF_nfe_2 = 0;
            }
        }else{
            AF_nfe_2 = AF_nfe->v1.f;
        }

        out_file_stream <<  AC_nfe_2 <<"\t";
        out_file_stream <<  AN_nfe_2 <<"\t";
        out_file_stream <<  AF_nfe_2 << "\t";

        //SOUTH ASIAN
        bcf_info_t* AN_sas  = bcf_get_info(hdr, rec, "AN_sas");
        bcf_info_t* AC_sas  = bcf_get_info(hdr, rec, "AC_sas");
        bcf_info_t* AF_sas  = bcf_get_info(hdr, rec, "AF_sas");
        int AC_sas_2 = AC_sas->v1.i; // allele count will always be an integer
        int AN_sas_2 = AN_sas->v1.i; // allele nufber will always be an interger

        float AF_sas_2;
        if(AF_sas == NULL) {
            if(AN_sas_2 != 0){
                AF_sas_2 = AC_sas_2 / AN_sas_2;
            }
            else{
                AF_sas_2 = 0;
            }
        }else{
            AF_sas_2 = AF_sas->v1.f;
        }

        out_file_stream << AC_sas_2 << "\t";
        out_file_stream << AN_sas_2 << "\t";
        out_file_stream << AF_sas_2 << "\t";


        bcf_info_t* cadd_phred  = bcf_get_info(hdr, rec, "cadd_phred");
        if(cadd_phred != NULL){
            int cadd_phred_2 = cadd_phred->v1.f; // allele count will always be an integer
            out_file_stream <<  cadd_phred_2 <<"\t";
        }else{
            out_file_stream << "NULL" <<"\t";
        }

        bcf_info_t* revel_max  = bcf_get_info(hdr, rec, "revel_max");
        if(revel_max != NULL){
            int revel_max_2 = revel_max->v1.f; // allele count will always be an integer
            out_file_stream  << revel_max_2 <<"\t";
        }else{
            out_file_stream << "NULL" <<"\t";
        }

        bcf_info_t* spliceai_ds_max = bcf_get_info(hdr, rec, "spliceai_ds_max");
        if(spliceai_ds_max != NULL){
            int spliceai_ds_max_2 = spliceai_ds_max->v1.f; // allele count will always be an integer
            out_file_stream << spliceai_ds_max_2 <<"\t";
        }else{
            out_file_stream << "NULL" <<"\t";
        }

        bcf_info_t* pangolin_largest_ds  = bcf_get_info(hdr, rec, "pangolin_largest_ds");
        if(pangolin_largest_ds != NULL){
            int pangolin_largest_ds_2 = pangolin_largest_ds->v1.f; // allele count will always be an integer
            out_file_stream << pangolin_largest_ds_2 <<"\t";
        }else{
            out_file_stream << "NULL" <<"\t";
        }

        bcf_info_t* phylop  = bcf_get_info(hdr, rec, "phylop");
        if(phylop != NULL){
            int phylop_2 = phylop->v1.f; // allele count will always be an integer
            out_file_stream << phylop_2 <<"\t";
        }else{
            out_file_stream << "NULL" <<"\t";
        }

        bcf_info_t* sift_max  = bcf_get_info(hdr, rec, "sift_max");
        if(sift_max != NULL){
            int sift_max_2 = sift_max->v1.f; // allele count will always be an integer
            out_file_stream <<  sift_max_2 <<"\t";
        }else{
            out_file_stream << "NULL" <<"\t";
        }

        bcf_info_t* polyphen_max  = bcf_get_info(hdr, rec, "polyphen_max");
        if(polyphen_max != NULL){
            int polyphen_max_2 = polyphen_max->v1.f; // allele count will always be an integer
            out_file_stream << polyphen_max_2 <<"\t";
        }else{
            out_file_stream << "NULL" <<"\t";
        }

        out_file_stream << vep_val;

        out_file_stream << std::endl;

    }

    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    bcf_destroy(rec);
}

