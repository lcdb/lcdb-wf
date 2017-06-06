# Demo wrapper

This wrapper demonstrates current best-practices.

The target audience of the wrapper's README should be yourself six months from
now, under a tight deadline, frantically looking for that rule you wrote so you
can copy/paste into a custom Snakefile.

Examples should come first. There should be at least a minimal example and
a reasonably complicated example. To be complete you can add links to docs,
a brief description of the tool, and example output.

This demo wrapper simply copies input files to output files.

## Examples

Minimal usage:

```python
rule demo:
    input: 'a.txt'
    output: 'b.txt'
    wrapper:
        'file://path/to/wrapper'
```

"paired-end" usage:

```python
rule demo:
    input:
        R1='a1.txt',
        R2='a2.txt'
    output:
        R1='b1.txt',
        R2='b2.txt'
    wrapper:
        'file://path/to/wrapper'
```

## Input

Input file formats for this wrapper can be anything.

### Single-end mode:

Expects a single unnamed input file.

### Paired-end mode:

Expects two input files with keys `R1` and `R2`.

## Output

Output files are simply copies of input.

### Single-end mode:

Expects a single unnamed output file

### Paired-end mode:

Expects two output files with keys `R1` and `R2`.

## Threads
Does not use threads

## Params
Does not use params
