import csv

#config: "config/config.json"
workdir: "/blue/kgraim/leslie.smith1/Repositories/PhyloFrameGNN/"
log_dir = "logs"

# create list of chromosome numbers we will be using (used as wildcards)
init_chrom = list(range(1,23))
CHROMOSOMES=list(map(str, init_chrom))
#CHROMOSOMES = list(range(1,23))
CHROMOSOMES.append("X")
CHROMOSOMES.append("Y")

# here_path should be the path to "project home"
here_path = "/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/PhyloFrameGNN/"
rule all:
    input:
        #expand(here_path + "processed-data/gnomADv4.1_EAF/Joint/gnomad.joint.v4.1.sites.chr{chrom_files}.eaf.tsv", chrom_files=CHROMOSOMES)
        expand(here_path + "processed-data/gnomADv4.1_EAF/Exomes/gnomad.exomes.v4.1.sites.chr{chrom_files}.eaf.tsv", chrom_files=CHROMOSOMES)

# download vcf files for each chromosome from gnomAD
rule download_chromosomes:
    #source script: /orange/kgraim/data/gnomAD/download_script.sh
    params:
        #url = "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{chrom_files}.vcf.bgz",
        url = "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr{chrom_files}.vcf.bgz",
        out_path = here_path + "raw-data/gnomADv4.1/Exomes_Raw/"
       
    log:
        stdout=log_dir+ "/" + "chr{chrom_files}" + "_" + "Download.stdout", 
        stderr=log_dir + "/" + "chr{chrom_files}" + "_" + "Download.stderr"
    output:
        out_file = here_path + "raw-data/gnomADv4.1/Exomes_Raw/gnomad.exomes.v4.1.sites.chr{chrom_files}.vcf.bgz"
        
    shell:
        """
        wget {params.url} -P {params.out_path}
        """

# parse ancestry frequencies for each chromosome
rule parse_chromosomes: 
    envmodules:
        "samtools",
        "cmake",
        "gcc/9.3.0"
    #params:
     #   out_path = here_path + "raw-data/gnomADv4.1/Exomes_Parsed/"
    log:
        stdout=log_dir+ "/" + "chr{chrom_files}" + "_" + "Parse.stdout", 
        stderr=log_dir + "/" + "chr{chrom_files}" + "_" + "Parse.stderr"
    input:
        download_file = here_path + "raw-data/gnomADv4.1/Exomes_Raw/gnomad.exomes.v4.1.sites.chr{chrom_files}.vcf.bgz"

    output:
        parsed = here_path + "raw-data/gnomADv4.1/Exomes_Parsed/gnomad.exomes.v4.1.sites.chr{chrom_files}.parsed.vcf" 
    shell:
        """
        /home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/PhyloFrame/vcf_parsing/build/read_freq -i {input.download_file} -o {output.parsed}
        """

# calculate enhanced allele freuencies for each chromosome
rule calculate_EAF:
    envmodules:
        "R"
       
    log:
        stdout=log_dir+ "/" + "chr{chrom_files}" + "_" + "EAF.stdout", 
        stderr=log_dir + "/" + "chr{chrom_files}" + "_" + "EAF.stderr"
    input:
        parsed = here_path + "raw-data/gnomADv4.1/Exomes_Parsed/gnomad.exomes.v4.1.sites.chr{chrom_files}.parsed.vcf" 
    output:
        eaf = here_path + "processed-data/gnomADv4.1_EAF/Exomes/gnomad.exomes.v4.1.sites.chr{chrom_files}.eaf.tsv"
    shell:
        """
        Rscript /blue/kgraim/leslie.smith1/Repositories/PhyloFrameGNN/code/gnomadV4_EAF_Calculation.R -i {input.parsed} -o {output.eaf}
        """
