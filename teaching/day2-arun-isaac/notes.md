You've all run multi-step bioinfo pipelines. What goes wrong? You rerun things unnecessarily. You forget which parameters you used. Steps fail and you can't tell which one. That's what CWL + ravanan solves.

# Part 0: Install the tools

We use Guix to drop into a shell with all the necessary tools.
```
guix shell -m manifest.scm
```

# Part 1: First CWL workflow

For our first CWL workflow, let's create a really simple one that just prints a given message. Put the following in a file `hello.scm`.
```
(define print
  (command #:inputs (message #:type string)
           #:run "echo" message
           #:outputs (output-message #:type stdout))

(workflow ((message #:type string))
  (print #:message message))
```
Compile it to CWL.
```
ccwl compile hello.scm -o hello.cwl
```
Set up an `hello-inputs.yaml` file.
```
message: foo
```
Then, run it.
```
ravanan --guix-manifest=ravanan-manifest.scm --store=store hello.cwl hello-inputs.yaml
```
It's also informative visualizing the graphical structure of the workflow. Try:
```
ccwl compile hello.scm -t dot | dot -Tpng -o hello.png
```

Now, let's try a slightly more complex two-step workflow—one that prints the sequence labels in a compressed FASTA file. Put the following in `fasta-sequences.scm`.
```
(define uncompress
  (command #:inputs (compressed-fasta #:type File)
           #:run "gunzip" "--decompress" "--to-stdout" compressed-fasta
           #:outputs (fasta #:type stdout)))

(define sequence-labels
  (command #:inputs (fasta #:type File)
           #:run "grep" "^>" fasta
           #:outputs (labels #:type stdout)))

(workflow ((compressed-fasta #:type File))
  (pipe (uncompress #:compressed-fasta compressed-fasta)
        (sequence-labels #:fasta fasta)))
```
The first step—the uncompress step—uncompresses the FASTA file. The second step—the sequence-labels step—prints out the sequence labels. The final workflow step connects the two steps to execute one after another.

Then we download a compressed FASTA file.
```
wget https://github.com/pangenome/pggb/raw/refs/heads/master/data/HLA/DRB1-3123.fa.gz
```
And, we set up an `fasta-sequences-inputs.yaml` file.
```
compressed-fasta:
  class: File
  path: DRB1-3123.fa.gz
```
Finally, we compile and run like last time.
```
ccwl compile fasta-sequences.scm -o fasta-sequences.cwl
ravanan --guix-manifest=ravanan-manifest.scm --store=store fasta-sequences.cwl fasta-sequences-inputs.yaml
```
Examine the output and see if it makes sense. Look at the visualization of the workflow structure too.
```
ccwl compile fasta-sequences.scm -t dot | dot -Tpng -o fasta-sequences.png
```
Run the workflow again without changing anything. What happens? Both steps should complete instantly. ravanan recognizes that nothing changed and reuses the cached results.

## Further reading

* **ccwl manual:** https://ccwl.systemreboot.net/manual/dev/en/
* **ravanan repo:** https://github.com/arunisaac/ravanan/

# Part 2: pggb.cwl
In this part, we will run pggb.cwl, a pangenome building workflow. pggb.cwl is a CWL port of pggb, the pangenome graph builder script that you will encounter later in this course.

First, clone the repo.
```
git clone https://github.com/arunisaac/pggb.cwl
cd pggb.cwl
```
Compile the workflow to CWL.
```
ccwl compile pggb.scm -o pggb.cwl
```
Download sequences from an example dataset in the pggb repo, and uncompress it.
```
wget https://github.com/pangenome/pggb/raw/refs/heads/master/data/HLA/DRB1-3123.fa.gz
gunzip DRB1-3123.fa.gz
```
Set up an `inputs.yaml` file.
```
sequences:
  class: File
  path: DRB1-3123.fa
number-of-haplotypes: 12
threads: 8
```
Then, run the workflow.
```
ravanan --guix-channels=channels.scm --store=store pggb.cwl inputs.yaml
```
Congratulations, you have just built your first pangenome! Look at the various files listed in the output. The pangenome visualizations are particularly pretty to look at. You will learn to interpret them later in the course.

## Further exploration
Look in `pggb.scm` (towards the end in the final workflow section) to find other supported parameters. Experiment with varying parameters. For example, try setting `map-percent-identity` to 95. You will learn what each of these parameters do later in the course.

Notice how only steps affected by the parameter changes are rerun.
