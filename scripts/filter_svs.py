import pysam
path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/"
vcf_in = path + 'structural/CBFH00195_sniffles.vcf'
vcf_out = path + 'processed/CBFH00195_sniffles_filteredremove.vcf'
txt_out = path + "processed/CBFH00195_sniffles_filteredremove.list"
threshold = 2  # Set your difference threshold
meancov = 30




with pysam.VariantFile(vcf_in) as vcf, open(vcf_out, "w") as out_vcf, open(txt_out, "w") as txt:
    out_vcf.write(str(vcf.header))  # Write header
    for record in vcf:
        coverage_str = record.info.get("COVERAGE")
        if coverage_str:
            coverage = list(map(int, coverage_str))
            if coverage[0] > 0 and coverage[4] > 0:
                if coverage[0] > 2*meancov or coverage[4]> 2*meancov:
                    cov_diff = coverage[0] / coverage[4] 
                    cov_diff2 = coverage[4] / coverage[0] 
                    if cov_diff > threshold or cov_diff2 > threshold:
                        out_vcf.write(str(record))
                        txt.write(str(record.id)+"\n")
            else: 
                out_vcf.write(str(record))
                txt.write(str(record.id)+"\n")
        type = record.info.get("SVTYPE")
        support = record.info.get("SUPPORT")
        if type == "BND" or support < meancov/4:
           out_vcf.write(str(record))
           txt.write(str(record.id)+"\n")           

        for sample in record.samples:
            genotype = record.samples[sample]['GT']
            if genotype == (0, 0):
                out_vcf.write(str(record))
                txt.write(str(record.id)+"\n")

### this script must be followed by
## grep -F -v -f CBFH00251_sniffles_filteredremove.list ../structural/CBFH00251_sniffles.vcf > CBFH00251_sniffles_clean.vcf

