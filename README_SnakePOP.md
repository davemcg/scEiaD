README-SnakePOP
================

SnakePOP is Snakefile for running batch integrations on large scRNA-seq.
This snakefile was written to allow different sets of integration
parameters - batch correction algorithmn, normalization method, etc - to
be added and run effciently. It allows for creating wildcards specific
to subset of data or a specific integration algorithmn

## Non-snakemake infrastructure

There’s a fair bit of code outside the actual snakemake rules to allow
for flexibility within SnakePOP, which is described below. While this
seems like a lot, it runs fast - almost no noticable change in
performance.

### Integration configuration

The most important input is the `integration-config`. I’ll reference the
file `integration-config-scVI_smalll.json` below.

    {
        "global_default" :{
            "covariate" : ["batch"]
        },
        "TabulaDroplet":{
            "default": {
                    "method" : ["scVI", "ldvae"],
                    "transform" : ["counts"],
                "n_features" : [5000],
                "dims" : [10],
                "knn" : [5],
                "dist": [0.1],
                "neighbors": [500]
            }
        },
    ...

The file is written in json format. All terminal entries(covariate,
n\_features etc.) represent integration parameters. The names of the
entries are used as wildcards within SnakePOP, and the values are the
values filled in for each wildcard. The top most tag, `global_default` ,
sets wildcards that will be run across *all* different combinations of
parameters. Tags below(but at the same level as) `global_default`
represent different data partitions. In the example above,
`TabulaDroplet` is a data partition the corresponds only to droplet
based data within scEiaD. Within each parition tag, there must be a
`default` entry. Wihtin the `default` tag, you specify the wildcards and
values that will be run specific to the partition. Tags below(but at the
same level as) `default` can be added to specify wildcards+ values
specific to one of the methods specified in `default.method` ie.

``` 
    "TabulaDroplet":{
        "default": {
                "method" : ["scVI", "ldvae"],
                "transform" : ["counts"],
            "n_features" : [5000],
            "dims" : [10],
        },
        "ldave" :{
            "n_features" : [1000],
            "epochs" : [10,20]
        }
    },
```

These method-specific values will supercede those in default. Any
parameter not specified within the method-specfic tag will inherit from
the `default` tag. Using this format, you can write out number of
partition specific and method specific combos. See
`integration-config-all.json` for the fully fleshed out integrations run
for scEiaD

### Rule input/output abstraction

In order to make it easy to add/remove wildcard parameters within the
Snakemake workflow, SnakePOP automatically generates input and output
strings containing specific wildcards for each rule. These strings are
generated with the function `make_wc_string`, and has three key
arguments:

  - pfx : a string for the location prefix for where to store the output
    file
  - wildcards: a list of wildcards to use. These *must* match the
    parameters specified in `integration-config`
  - sfx : a string for the location suffix for the output file

example(for a rule within no input)

    exp_pfx = 'data/example/
    exp_wc = ['method', 'dims']
    exp_sfx = '.txt'
    rule example:
        output: make_wc_string(exp_pfx, exp_wc, exp_sfx)

The same `make_wc_string` call can be used as the input for downstream
rules. As new wildcards are added in downstream rules, the inputs for
make\_wc\_string can be changed accordingly. The best way to do this is
to append to a previously specfied list of wildcards. This helps trace
out the flow of wildcards throughout the pipeline, and ensures that
wildcard values specified at the end of the pipeline are passed back
properly.

    downstream_pfx = 'data/integrated_data/'
    downstream_wc = exp_wc +  ['Epochs', 'layers]
    downstream_sfx = '.Rdata'
    rule downstream:
        input: make_wc_string(exp_pfx, exp_wc, exp_sfx) 
        output: make_wc_string(downstream_pfx, downstream_wc, downstream_sfx)

It’s useful to trace out the order in which wildcards are inherited and
added throught the workflow, which looks like this

    seu_obj_param|
                 |-> intg_obj_param == phate_obj_param
                                   |-> umap_obj_param| 
                                                     |->extr_umap_param ==  perf_met_param == scIB_param 
                                   |-> tsne_obj_param == extr_tsne_param
                                   |-> clst_obj_param == extr_clst_param 

In the final rule, `merge_stats`, values for the wildcards are filled in
using the function `parse_integration_config`. This function has three
inputs - file: the path to the `integration-config` file. Right now, its
specified in `config.yaml` - output\_string: a string generated by
`make_wc_string` the will be filled in with wildcards from -
which\_partitions: a set of paritions specified within
`integration-config` that will be run.

parse\_integration\_config will return a list where each entry
corresponds to a specific combination of parameters within the
`integration-config`.

## Rule level run parameter export

To make it easy to pass specific wildcard values associated with each
run of a rule to an Rscript, we can export the input, output, and
wildcards for a specific run of a rule as a yaml, using the function
`export_rule_info`. These yamls are written to a folder , `RSON_temp`
which is specified in the config.yaml.

    input:
        ...
    output:
        ...
    params:
        rson = lambda wildcards, input, output: export_rule_info(input = dict(input), output=dict(output), wildcards = dict(wildcards))

`export_rule_info` uses kwargs, so any dict can be passed as a named
argument

## Rules in main worflow

### inputs to scPOP

The two biggest inputs are a genes x cells sparseMatrix Rdata file, and
cell level meta-data file (`cell_info`), where each
barcode-sample\_accession id corresponds to a columns in the
sparseMatrix.

### rule label\_known\_cells\_with\_type

  - reads in metadata from studies that have labelled data, and adds
    labels to `cell_info`.

### rule make\_seurat\_obj

  - calls `src/build_seurat_obj_classic.R`. splits up data based on
    `parition` wildcards, and if required normalizes data. Most the
    important code is sourced from `src/make_seurat_obj_functions.R` -
    this handles much of the normalization of differnt partitions/
    technology types.

### rule integrate00

  - runs the different integration algorithmns. calls
    `src/merge_methods.R`, which configures conda environments properly
    to run the different tools
