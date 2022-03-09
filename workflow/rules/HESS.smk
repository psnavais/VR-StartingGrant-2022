rule format_sumstats_HESS:
        'Format summary statistics for HESS heritability estimation.'
        input:
                '/mnt/work2/pol/other_pheno/miscarriage/results/delivery/{sample}/{pheno}.allchr.txt.gz'
        output:
                '/mnt/work2/pol/other_pheno/miscarriage/HESS/aux/sumstats/{pheno}-{sample}.txt'
        run:
                d= pd.read_csv(input[0], sep= ' ', header= 0, usecols= ['CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 'N', 'BETA', 'SE', 'A1FREQ'])[['CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 'N', 'BETA', 'SE', 'A1FREQ']]
		d['MAF']= np.where(d.A1FREQ> 0.5, 1 -d.A1FREQ, d.A1FREQ )
		d= d.loc[d.MAF> 0.05, :]
		d['SNP']= np.where(d.ALLELE1> d.ALLELE0, d.CHROM.apply(str) + ':' + d.GENPOS.apply(str) + ':' + d.ALLELE0 + ':' + d.ALLELE1, d.CHROM.apply(str) + ':' + d.GENPOS.apply(str) + ':' + d.ALLELE1 + ':' + d.ALLELE0)
                d.drop_duplicates(['CHROM', 'GENPOS'], inplace= True)
                d['Z']= d.BETA / d.SE
                d.columns= ['CHR', 'BP', 'A2', 'A1', 'N', 'BETA', 'SE', 'A1FREQ', 'MAF', 'SNP', 'Z']
                d.to_csv(output[0], sep= '\t', header= True, index= False, columns= ['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z', 'N'])

rule HESS_step1:
        'Run HESS step1.'
        input:
                '/mnt/work2/pol/other_pheno/miscarriage/HESS/aux/sumstats/{pheno}-{sample}.txt',
                '/mnt/work2/pol/refdata/ld_indep_regions.txt',
                multiext('/mnt/work2/pol/refdata/1kg_eur_1pct/1kg_eur_1pct_chr{autoCHR}', '.bim', '.bed', '.fam')
        output:
                temp(multiext('/mnt/work2/pol/other_pheno/miscarriage/HESS/step1/{pheno}-{sample}_chr{autoCHR}', '.info.gz', '.eig.gz', '.prjsq.gz', '.log'))
        params:
                '/mnt/work2/pol/refdata/1kg_eur_1pct/1kg_eur_1pct_chr{autoCHR}',
                '/mnt/work2/pol/other_pheno/miscarriage/HESS/step1/{pheno}-{sample}'
        conda:
                './HESS.yaml'
        shell:
                '''
                python2 ~/soft/hess-0.5.3-beta/hess.py \
                --local-hsqg {input[0]} \
                --chrom {wildcards.autoCHR} \
                --bfile {params[0]} \
                --partition {input[1]} \
                --out {params[1]}
                '''
rule HESS_step2:
        'Run HESS step2.'
        input:
                expand('/mnt/work2/pol/other_pheno/miscarriage/HESS/step1/{{pheno}}-{{sample}}_chr{autoCHR}.info.gz', autoCHR= autosomal_CHR),
                expand('/mnt/work2/pol/other_pheno/miscarriage/HESS/step1/{{pheno}}-{{sample}}_chr{autoCHR}.eig.gz', autoCHR= autosomal_CHR),
                expand('/mnt/work2/pol/other_pheno/miscarriage/HESS/step1/{{pheno}}-{{sample}}_chr{autoCHR}.prjsq.gz', autoCHR= autosomal_CHR),
                expand('/mnt/work2/pol/other_pheno/miscarriage/HESS/step1/{{pheno}}-{{sample}}_chr{autoCHR}.log', autoCHR= autosomal_CHR)
        output:
                multiext('/mnt/work2/pol/other_pheno/miscarriage/HESS/step2/{pheno}-{sample}', '.log', '.txt')
        params:
                '/mnt/work2/pol/other_pheno/miscarriage/HESS/step1/{pheno}-{sample}',
                '/mnt/work2/pol/other_pheno/miscarriage/HESS/step2/{pheno}-{sample}'
        threads: 10
        conda:
                './HESS.yaml'
        shell:
                '''
                python2 ~/soft/hess-0.5.3-beta/hess.py --prefix {params[0]} --out {params[1]}
                '''

rule check_HESS:
        'Check HESS is created.'
        input:
                expand('/mnt/work2/pol/other_pheno/miscarriage/HESS/step2/{pheno}-{sample}.{ext}', ext= ['txt'], pheno= pheno_nms, sample= sample_nms)
        output:
                '/mnt/work2/pol/other_pheno/miscarriage/HESS/checks/h2_estimated.txt'
        shell:
                'touch {output[0]}'

