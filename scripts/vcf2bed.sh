grep -v "#" CBFH00251_sniffles_clean.vcf | awk 'BEGIN { OFS="\t" }
      	!/^#/ {
      	split($8, info, ";");
                 end = ""; svtype = "";
                 for (i in info) {
                     if (info[i] ~ /^END=/) {
                         end = substr(info[i], 5);
                     } 
                     if (info[i] ~ /^SVTYPE=/) {
                         svtype = substr(info[i], 8);
                      }
                  }
                  print $1, $2 - 1, end, $3, $10
              }' - > CBFH00251_sniffles_clean.bed


grep -v "#" CBFH00195_sniffles_clean.vcf | awk 'BEGIN { OFS="\t" }
      	!/^#/ {
      	split($8, info, ";");
                 end = ""; svtype = "";
                 for (i in info) {
                     if (info[i] ~ /^END=/) {
                         end = substr(info[i], 5);
                     } 
                     if (info[i] ~ /^SVTYPE=/) {
                         svtype = substr(info[i], 8);
                      }
                  }
                  print $1, $2 - 1, end, $3, $10
              }' - > CBFH00195_sniffles_clean.bed

