#!/usr/bin/env bash

# 1. Set variables to make the script works from any directory
BASH_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
SM_DIR=$BASH_DIR/../..
cd $SM_DIR

OUTDIR=out/wget_interesting_samples_from_GSE60104/
mkdir -p $OUTDIR
cd $OUTDIR

# 2. Query produced by Cliget Firefox add-on after selecting the 4 desired samples from custom download here:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60104
curl --header 'Host: www.ncbi.nlm.nih.gov' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:81.0) Gecko/20100101 Firefox/81.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60104' --header 'Content-Type: application/x-www-form-urlencoded' --header 'Origin: https://www.ncbi.nlm.nih.gov' --header 'DNT: 1' --cookie 'ncbi_sid=0C422A67B4A055C1_0000SID; WebCubbyUser=HQR07QIUYXRKHJQ9RHJT8Z9RCD33NDSJ%3Blogged-in%3Dtrue%3Bmy-name%3Dguillaumecharbonnier%3Bpersistent%3Dtrue%400C422A67B4A055C1_0000SID; entrezSort=pubmed:none; pmc.article.report=; books.article.report=; MyNcbiSigninPreferences=O25jYmlscyY%3D; ncbi_pinger=N4IgDgTgpgbg+mAFgSwCYgFwgAwGEAsATAJwBC++AYgBwAipuArNQIzbsedcsgA0IAYwA2yAQGsAdlAAeAF0ygeWANoB7AEYArKANkACABIAVALIAZAJISwAV1kBRIVAC2UCbIC6vATYDOs1WdaIQBBaABDPhAAZmJMGIA2FkZo/Cj8bHjohNTiAHZ0tKxhUUkZeX58RniCyoT4gDNwoV8odIKsfDzo9Op4wm7qyrisKMZM0f5GJUE/AKChAGUbdWdkCpBGQnjaVQB3CSFVcPQAX34bQ+PUKTkFGJ6sJpa2/lT42QgbV5i+yZiRolkvg+pUJolotNtpUiiBLkcTrcNlV4qDNvUsCRapsOiA2HF+El4gkiCBzoJAs5VBIkfdtlgAOZQVRRR4gdKNGxCIRjeJRQjgvBEMgUGj0JisLhSzg8fiEGZM1QYAByAHllfYMA0uUIVer7Pz6SA9iaAHQSATqZDmoTOc3IRCmhmqGD8wEsYgC1ngtjUTJvGYer1vI1sLbpNnYU3RKM9KZs4jMMawz09U6nIA=; WebEnv=1Fiowb%400C422A67B4A055C1_0000SID' --header 'Upgrade-Insecure-Requests: 1' --request POST --data-urlencode 'acc=GSE60104' --data-urlencode 'format=file' --data-urlencode 'id=1,19,27,59' 'https://www.ncbi.nlm.nih.gov/geo/download/' --output 'GSE60104_RAW.tar'

tar xf GSE60104_RAW.tar
rm -f GSE60104_RAW.tar
