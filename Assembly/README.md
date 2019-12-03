## Parameter for FALCON genome assembly
### FALCON only assembly for 75-fold NC358:
>pa_HPCdaligner_option =  -k14 -e0.75 -s100 -l3000 -h240 -w8 -H14154
ovlp_HPCdaligner_option =  -k24 -e.95 -s100 -l1000 -h600 -H20149
pa_DBsplit_option = -x500 -s400
ovlp_DBsplit_option = -x500 -s400
falcon_sense_option = --min_idt 0.70 --min_cov 2 --max_n_read 200
overlap_filtering_setting = --max_diff 40 --max_cov 80 --min_cov 2 --min_len 500

### Canu only assembly for 75-fold NC358:
>ovlMerThreshold=500; genome_size=2272400000; input_type=-pacbio-raw

## FALCON-Canu hybrid assembly for 60-fold NC358:
### FALCON
>pa_HPCdaligner_option =  -k14 -e0.75 -s100 -l3000 -h240 -w8 -H10310
pa_DBsplit_option = -x500 -s400  
falcon_sense_option = --min_idt 0.70 --min_cov 2 --max_n_read 200
### CANU
>ovlMerThreshold=500; genome_size=2272400000; input_type=-pacbio-corrected

## FALCON-Canu hybrid assembly for 50, 40, 30, 20-fold NC358:
### FALCON
>pa_HPCdaligner_option =  -k14 -e0.75 -s100 -l3000 -h240 -w8 -H3000
pa_DBsplit_option = -x500 -s400  
falcon_sense_option = --min_idt 0.70 --min_cov 2 --max_n_read 200
### CANU
>ovlMerThreshold=500; genome_size=2272400000; input_type=-pacbio-corrected

## FALCON-Canu hybrid assembly for 50-fold NC358 with shift distribution:
### FALCON
>pa_HPCdaligner_option =  -k14 -e0.75 -s100 -l3000 -h240 -w8 -H3000
pa_DBsplit_option = -x500 -s400  
falcon_sense_option = --min_idt 0.70 --min_cov 2 --max_n_read 200
### CANU
>ovlMerThreshold=500; genome_size=2272400000; input_type=-pacbio-corrected


## FALCON only assembly for 68-fold B73:
>pa_HPCdaligner_option =  -k14 -e0.75 -s100 -l3000 -h240 -w8 -H9898
ovlp_HPCdaligner_option =  -k29 -e.95 -s100 -l4800 -h600 -H12360
pa_DBsplit_option = -x500 -s400  
ovlp_DBsplit_option = -x500 -s400
falcon_sense_option = --min_idt 0.70 --min_cov 2 --max_n_read 200
overlap_filtering_setting = --max_diff 40 --max_cov 80 --min_cov 2 --min_len 500

## Canu only assembly for 68-fold B73:
>ovlMerThreshold=500; genome_size=2500000000; input_type=-pacbio-raw

## FALCON-Canu hybrid assembly for 68-fold B73:
### FALCON
>pa_HPCdaligner_option =  -k14 -e0.75 -s100 -l3000 -h240 -w8 -H9898
pa_DBsplit_option = -x500 -s400  
falcon_sense_option = --min_idt 0.70 --min_cov 2 --max_n_read 200
### CANU
>ovlMerThreshold=500; genome_size=2272400000; input_type=-pacbio-corrected

