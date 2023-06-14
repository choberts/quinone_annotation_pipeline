from glob import glob

configfile: "config.yml"

DIR = "data/"+config["dataset"]
genomes_dir = {genome: os.path.basename(genome) for genome in glob(DIR+"/GC*")}


rule all: 
    input:
        expand("results/{dataset}/{dataset}_annot_res.tsv", dataset=config["dataset"])
        #expand("results/{dataset}/{dataset}_hmmscan_agg.tsv", dataset=config["dataset"])
        #expand("results/{dataset}/hmmscan_output/{id}_hmmscan.tsv", dataset=config["dataset"], id=genomes_dir.values()) 


rule unzip:
     input:
         "data/{dataset}/{id}/{id}_protein.faa.gz"
     output:
         temporary("data/{dataset}/{id}/{id}_protein.faa")
     shell:
         "gzip -d {input} --keep"


rule run_hmmscan:
    params:
        profile = config['hmm_profiles']
    input:
        "data/{dataset}/{id}/{id}_protein.faa"
    output: 
        "results/{dataset}/hmmscan_output/{id}_hmmscan.tsv"
    shell:
        "hmmscan --domtblout {output} {params.profile} {input}"


rule aggregate_hmmscan:
    params: 
        id = list(genomes_dir.values()),
        dataset = config['dataset'] 
    input:
        expand("results/{{dataset}}/hmmscan_output/{id}_hmmscan.tsv", id=genomes_dir.values())
    output:
        "results/{dataset}/{dataset}_hmmscan_agg.tsv"
    shell:
        "for el in {params.id}; do G=$el; cat results/{params.dataset}/hmmscan_output/$G\_hmmscan.tsv | tail -n +4 | head -n -10 | awk -v GCF=$G '/1/ {{print GCF, $0}}'; done > {output}" 


rule generate_tables:
    params:
        tax = config["taxonomy_path"]
    input:
        hmmscan_res = "results/{dataset}/{dataset}_hmmscan_agg.tsv",
        script = "scripts/parse_hmmscan_agg.py"
    output:
        "results/{dataset}/{dataset}_annot_res.tsv"
    shell:
        "python {input.script} -F {input.hmmscan_res} -T {params.tax} -O {output} -i 0.005 -c 0.5"

