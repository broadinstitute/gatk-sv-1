version 1.0

import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow GeneratePESRMetrics {
    input {
        File vcf
        String prefix

        File mean_coverage_file
        File? pe_file
        File? sr_file

        File medianfile
        File ped_file
        File male_samples
        File female_samples

        Int records_per_shard

        Int? pe_inner_window
        Int? pe_outer_window
        Int? sr_window
        String? additional_gatk_args

        Float? java_mem_fraction

        String gatk_docker
        String linux_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_scatter
        RuntimeAttr? runtime_attr_agg_pesr
        RuntimeAttr? runtime_override_concat
    }

    call tasks.ScatterVcf {
        input:
            vcf=vcf,
            records_per_shard = records_per_shard,
            prefix = "~{prefix}.scatter_vcf",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_scatter
    }

    scatter ( i in range(length(ScatterVcf.shards)) ) {
        call AggregatePESREvidence {
            input:
                vcf = ScatterVcf.shards[i],
                output_prefix = "~{prefix}.aggregate.shard_~{i}",
                mean_coverage_file = mean_coverage_file,
                pe_file = pe_file,
                sr_file = sr_file,
                pe_inner_window = pe_inner_window,
                pe_outer_window = pe_outer_window,
                sr_window = sr_window,
                additional_args=additional_gatk_args,
                java_mem_fraction = java_mem_fraction,
                gatk_docker = gatk_docker,
                runtime_attr_override = runtime_attr_agg_pesr
        }
    }

    call MiniTasks.ConcatVcfs as ConcatCleanedVcfs {
        input:
            vcfs=AggregatePESREvidence.out,
            vcfs_idx=AggregatePESREvidence.out_index,
            allow_overlaps=true,
            outfile_prefix="~{cohort_name}.aggregate_pesr",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_override_concat
    }

    output {
        File stats = MergeStats.merged_stats
        File? stats_common = MergeStatsCommon.merged_stats
    }
}


task AggregatePESREvidence {
    input {
        File vcf
        String output_prefix

        File mean_coverage_file
        File pe_file
        File sr_file

        Int? pe_inner_window
        Int? pe_outer_window
        Int? sr_window
        String? additional_args

        Float? java_mem_fraction
        String gatk_docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        pe_file: {
                     localization_optional: true
                 }
        sr_file: {
                     localization_optional: true
                 }
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcf, "GB") * 2.5),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.vcf.gz"
        File out_index = "~{output_prefix}.vcf.gz.tbi"
    }
    command <<<
        set -euo pipefail

        function getJavaMem() {
        # get JVM memory in MiB by getting total memory from /proc/meminfo
        # and multiplying by java_mem_fraction
        cat /proc/meminfo \
            | awk -v MEM_FIELD="$1" '{
                    f[substr($1, 1, length($1)-1)] = $2
                } END {
                    printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
                }'
            }
        JVM_MAX_MEM=$(getJavaMem MemTotal)
        echo "JVM memory: $JVM_MAX_MEM"

        gatk --java-options "-Xmx${JVM_MAX_MEM}" AggregatePairedEndAndSplitReadEvidence \
            -V ~{vcf} \
            -O ~{output_prefix}.vcf.gz \
            --sample-coverage ~{mean_coverage_file} \
            ~{"--discordant-pairs-file " + pe_file} \
            ~{"--split-reads-file " + sr_file} \
            ~{"--pe-inner-window " + pe_inner_window} \
            ~{"--pe-outer-window " + pe_outer_window} \
            ~{"--sr-window " + sr_window} \
            ~{additional_args}
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: gatk_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExtractMetrics {
    input {
        File vcf
        String output_prefix

        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcf, "GB") * 20),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.tsv"
    }
    command <<<
        set -euo pipefail
        echo "name PE_log_pval PE_called_median PE_bg_median PE_bg_frac SR_posA_log_pval SR_posB_log_pval SR_sum_log_pval SR_posA_called_median SR_posB_called_median SR_sum_called_median SR_posA_bg_median SR_posB_bg_median SR_sum_bg_median SR_posA_bg_frac SR_posB_bg_frac SR_sum_bg_frac SR_posA_pos SR_posB_pos PESR_log_pval PESR_called_median PESR_bg_median PESR_bg_frac" \
            | tr ' ' '\t' \
            > ~{output_prefix}.tsv
        bcftools query ~{vcf} \
            -f '%ID\t%PEQ\t%PECS\t%SR1Q\t%SR1CS\t%SR2Q\t%SR2CS\n"
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}