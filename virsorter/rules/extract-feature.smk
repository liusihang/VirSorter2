rule gff_feature:
    input: f'{Tmpdir}/all.pdg.gff'
    output: f'{Tmpdir}/all.pdg.gff.ftr'
    log: f'log/{Tmpdir}/step2-extract-feature/extract-feature-from-gff-common.log'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        python {Scriptdir}/extract-feature-from-gff.py {Dbdir}/rbs/rbs-catetory.tsv {input} {output} &> {log} || {{ echo "See error details in {Wkdir}/{log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """

rule gff_feature_by_group:
    input:
        ftr=f'{Tmpdir}/all.pdg.gff.ftr',
        gff=f'{Tmpdir}/{{group}}/all.pdg.gff'
    output: f'{Tmpdir}/{{group}}/all.pdg.gff.ftr'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step2-extract-feature/extract-feature-from-gff-{wildcards.group}.log
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        if [ -s $Rbs_pdg_db ]; then
            python {Scriptdir}/extract-feature-from-gff.py {Dbdir}/rbs/rbs-catetory.tsv {input.gff} {output} &> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            (cd {Tmpdir}/{wildcards.group} && ln -sf ../all.pdg.gff.ftr)
        fi
        """

##############pyhmmsearch##############################
rule hmmsearch:
    input:
        faa = f'{Tmpdir}/all.pdg.faa',
        faa_NCLDV = f'{Tmpdir}/NCLDV/all.pdg.faa'
    output: 
        tax = f'{Tmpdir}/all.pdg.hmm.tax',
        tax_NCLDV =f'{Tmpdir}/NCLDV/all.pdg.hmm.tax'
    params: 
        log = f'{Tmpdir}/hmmsearch.log' 
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        set -e # 如果任何命令失败，则立即停止执行
        python {Scriptdir}/pyhmmersearch.py --faa_file {input.faa} --db_dir {Dbdir} --output_file {output.tax} --evalue 1e-10 --bitscore 30 &> {params.log}
        python {Scriptdir}/selectbesthit.py {output.tax}
        python {Scriptdir}/selectbesthit.py {output.tax_NCLDV}
        """

######################Added feature for pyhmmersearch######################
rule hmm_sort_to_best_hit_taxon:
    input:
        faa = f'{Tmpdir}/{{group}}/all.pdg.faa',
        tax = f'{Tmpdir}/{{group}}/all.pdg.hmm.tax'
    output: 
        taxwhm = f'{Tmpdir}/{{group}}/all.pdg.hmm.taxwhm',
        ftr = f'{Tmpdir}/{{group}}/all.pdg.hmm.ftr'
    #log: f'log/{Tmpdir}/step2-extract-feature/extract-feature-from-hmmout-common.log'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {input.tax} > {output.ftr}

        # add hallmark info to .tax file for making affi-contigs.tab file
        python {Scriptdir}/add-hallmark-to-taxfile.py {input.tax} {output.taxwhm}
        # remove split faa files
        echo {Tmpdir}/all.pdg.faa.*.split | xargs rm -f
        """


rule hmmsearch_by_group:
    input:
        tax = f'{Tmpdir}/all.pdg.hmm.tax'
    output:
        output_tax = f'{Tmpdir}/{{group}}/all.pdg.hmm.tax'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        if [ -f {output.output_tax} ]; then
            echo "File {output.output_tax} already exists, skipping copy."
        else
            if [ -f {input.tax} ]; then
                cp {input.tax} {Tmpdir}/{wildcards.group}/ || {{ echo "Failed to copy {input.tax} to {Tmpdir}/{wildcards.group}/"; exit 1; }}
            else
                echo "Error: {input.tax} not found."
                exit 1
            fi
        fi
        """

localrules: hmm_features_by_group
rule hmm_features_by_group:
    input:
        ftr = f'{Tmpdir}/all.pdg.hmm.ftr',
        tax = f'{Tmpdir}/{{group}}/all.pdg.hmm.tax',
        faa = f'{Tmpdir}/{{group}}/all.pdg.faa'
    output: f'{Tmpdir}/{{group}}/all.pdg.hmm.ftr'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step2-extract-feature/merge-feature-{wildcards.group}.log
        Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
        Group_specific_hmmdb={Dbdir}/group/{wildcards.group}/customized.hmm
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db

        if [ -s $Hallmark_list_f ]; then
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {input.tax} --hallmark $Hallmark_list_f > {output} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            if [ {Prep_for_dramv} = True ]; then
                python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {input.tax}pfam --hallmark $Hallmark_list_f > {output}pfam 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi

        elif [ -s $Rbs_pdg_db ] || [ -s $Group_specific_hmmdb ]; then
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {input.tax} > {output} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }} 
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {input.tax}pfam > {output}pfam 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }} 
        else
            (cd {Tmpdir}/{wildcards.group} && ln -fs ../all.pdg.hmm.ftr)
            if [ {Prep_for_dramv} = True ]; then
                (cd {Tmpdir}/{wildcards.group} && ln -fs ../all.pdg.hmm.ftrpfam)
            fi
        fi
        """

localrules: merge_hmm_gff_features_by_group
rule merge_hmm_gff_features_by_group:
    input:
        gff_ftr = f'{Tmpdir}/{{group}}/all.pdg.gff.ftr',
        hmm_ftr = f'{Tmpdir}/{{group}}/all.pdg.hmm.ftr'
    output: 
        merged_ftr = f'{Tmpdir}/{{group}}/all.pdg.ftr'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step2-extract-feature/merge-feature-{wildcards.group}.log
        python {Scriptdir}/merge-hmm-gff-features.py {input.gff_ftr} {input.hmm_ftr} > {output.merged_ftr} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """
