.. _guide:

Guide to files
==============

The following is an annotated directory tree of the ``lcdb-wf`` repository to
help orient you. Hover over files for a tooltip description; click a file to
view the most recent version on GitHub.

Files in bold are the most important.

.. raw:: html

    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/balloon-css/0.2.4/balloon.min.css">

    <style>
    .dir {
        font-family: monospace;
        font-size: 1em;
    }
    .file {
        font-family: monospace;
        font-size: 0.8em;
    }
    .important {
        font-weight: bold;
    }
    .undoc {
        color: #888;
        }

    </style>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//README.md" style="text-decoration:none;"><span class="file undoc">README.md</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//references.snakefile" data-balloon=" The main workflow to generate references" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">references.snakefile</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//rnaseq.snakefile" data-balloon=" The main workflow for performing RNA-seq analysis" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">rnaseq.snakefile</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//4c.snakefile" data-balloon=" Worflow for 4C analysis using 4C-ker" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">4c.snakefile</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//requirements.txt" data-balloon=" Dependencies required for running lcdb-wf" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">requirements.txt</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//.travis.yml" style="text-decoration:none;"><span class="file undoc">.travis.yml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//.buildkite/" style="text-decoration:none;"><span class="dir undoc">.buildkite/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//.buildkite/pipeline.yml" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;pipeline.yml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//ci/" data-balloon=" Tools for managing the continuous integration tests" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="dir">ci/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//ci/build-docs.sh" data-balloon=" Builds documentation on travis-ci and automatically pushes to the gh-pages branch on github" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;build-docs.sh</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//ci/dependency_consistency.py" data-balloon=" Helper script for consistently updating dependencies of wrappers and requirements.txt." data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;dependency_consistency.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//ci/get-data.py" data-balloon=" Script for downloading example data." data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">&nbsp;&nbsp;&nbsp;get-data.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//ci/key.enc" data-balloon=" Encoded private key that allows pushing to github from travis-ci" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;key.enc</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//ci/travis-run.sh" data-balloon=" Runs tests on travis-ci" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;travis-run.sh</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//ci/travis-setup.sh" data-balloon=" Sets up environment on travis-ci" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;travis-setup.sh</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/" data-balloon=" This directory contains various configuration files used by the workflows" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="dir important">config/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/sampletable.tsv" data-balloon=" Table of sample metadata" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">&nbsp;&nbsp;&nbsp;sampletable.tsv</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/config.yml" data-balloon=" Main config file" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">&nbsp;&nbsp;&nbsp;config.yml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/clusterconfig.yaml" data-balloon=" Example cluster config file for running jobs on a cluster." data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">&nbsp;&nbsp;&nbsp;clusterconfig.yaml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/4c-sampletable.tsv" data-balloon=" Example sampletable for running a 4C analysis" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;4c-sampletable.tsv</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/envs/" data-balloon=" Conda environment definitions for per-rule environments that are not already a wrapper" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="dir">&nbsp;&nbsp;&nbsp;envs/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/envs/4c.yaml" data-balloon=" 4C workflow environment" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4c.yaml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/envs/annotationdbi.yaml" data-balloon=" Environment for building TSCs from AnnotationDbi tables" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;annotationdbi.yaml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/envs/references_env.yml" data-balloon=" Global environment for the references workflow" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;references_env.yml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/envs/R_rnaseq.yaml" data-balloon=" Environment for running RNA-seq differential expression" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;R_rnaseq.yaml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/multiqc_config.yaml" data-balloon=" Config file with additional settings for running MultiQC" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;multiqc_config.yaml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/test_4c_config.yaml" data-balloon=" Test config for 4C" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;test_4c_config.yaml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//config/test_config.yaml" data-balloon=" Test config for rnaseq" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;test_config.yaml</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//downstream/" style="text-decoration:none;"><span class="dir undoc">downstream/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//downstream/4c.R" data-balloon=" R script called by 4c.snakefile" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">&nbsp;&nbsp;&nbsp;4c.R</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//downstream/rnaseq.Rmd" data-balloon=" Rmd file called by rnaseq.snakefile" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">&nbsp;&nbsp;&nbsp;rnaseq.Rmd</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//downstream/rnaseq-requirements.txt" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;rnaseq-requirements.txt</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//include/" style="text-decoration:none;"><span class="dir undoc">include/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//include/adapters.fa" data-balloon=" Used in the cutadapt rules." data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;adapters.fa</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//include/WRAPPER_SLURM" data-balloon=" Wrapper script to submit jobs to a SLURM cluster" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file important">&nbsp;&nbsp;&nbsp;WRAPPER_SLURM</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/" data-balloon=" Directory of utilities used by the workflows" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="dir">lib/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/common.py" data-balloon=" The main module of utilities" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;common.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/fourc/" data-balloon=" 4C-specific utilities and scripts" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="dir">&nbsp;&nbsp;&nbsp;fourc/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/fourc/bedgraph_to_bigwig.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;bedgraph_to_bigwig.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/fourc/find-up-dn.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;find-up-dn.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/" data-balloon=" A package for post-processing references after they are downloaded." data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="dir important">&nbsp;&nbsp;&nbsp;postprocess/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/adapters.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;adapters.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/dicty.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dicty.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/dm6.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dm6.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/erccFisher.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;erccFisher.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/ercc.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ercc.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/hg19.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;hg19.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/hg38.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;hg38.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/__init__.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;__init__.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/merge.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;merge.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/phix.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;phix.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//lib/postprocess/sacCer3.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sacCer3.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//make_trackhub.py" style="text-decoration:none;"><span class="file undoc">make_trackhub.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/" style="text-decoration:none;"><span class="dir undoc">wrappers/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/.gitignore" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;.gitignore</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/LICENSE" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;LICENSE</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/README.md" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;README.md</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/" data-balloon=" Main test directory for wrappers" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="dir">&nbsp;&nbsp;&nbsp;test/</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/conftest.py" data-balloon=" Fixtures are imported here and used across py.test tests" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;conftest.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/raw_data_fixtures.py" data-balloon=" Fixtures for downloading example data" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="file">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;raw_data_fixtures.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_atropos.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_atropos.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_bowtie2.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_bowtie2.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_cutadapt.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_cutadapt.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_deeptools.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_deeptools.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_demo.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_demo.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_dupradar.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_dupradar.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_fastqc.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_fastqc.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_fastq_screen.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_fastq_screen.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_featurecounts.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_featurecounts.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_hisat2.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_hisat2.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_kallisto.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_kallisto.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_multiqc.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_multiqc.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_picard.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_picard.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_rseqc.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_rseqc.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_salmon.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_salmon.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/test_samtools.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test_samtools.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test_toy.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;test_toy.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/test/utils.py" style="text-decoration:none;"><span class="file undoc">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;utils.py</span></a></p>
    
    <p style="margin:0px;"><a href="https://github.com/lcdb/lcdb-wf/blob/master//wrappers/wrappers/" data-balloon=" Wrappers directory of snakemake wrappers used among the workflows" data-balloon-pos="right" data-balloon-length="xlarge" style="text-decoration:none;"><span class="dir important">&nbsp;&nbsp;&nbsp;wrappers/</span></a></p>
    


