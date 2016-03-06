# clj-tandem

A parser for X! Tandem XML formatted result files.

## Usage

Import from Clojars:

```clojure
[clj-tandem "0.1.1"]
```

Use in your namespace:

```clojure
(:require [clj-tandem.core :as xt])
```

Open a reader on a X! Tandem result file containing and call
'group-seq'. This will return a lazy list of zippers, one for each
model in the file, that can be used with the usual Clojure XML parsing
libraries. 'peptide-seq' called on a group zipper returns a collection
of maps containing information about the peptide identification.

```clojure
clj-tandem.core> (def xtf "/path/to/file.xml")
#'clj-tandem.core/xtf
clj-tandem.core> (with-open [r (reader tf)]
                   (->> (group-seq r)
                        (mapcat peptide-seq)
                        first))
{:protein-info {:expect "-34.9", :id "2693.1", :uid "10067", :label
"tr|I3LJX2|I3LJX2_PIG Uncharacterized protein OS=Sus scrofa PE=1 SV=1",
:sumI "5.21", :description "Uncharacterized protein OS=Sus scrofa PE=1 SV=1",
:accession "tr|I3LJX2|I3LJX2_PIG"}, :mh "1099.5491", :mods (), :pre "AGPK",
:expect 0.017, :protein-length 221, :missed_cleavages "1", :rank 1, :start 137,
:x_ions "0", :b_ions "3", :nextscore "353", :hyperscore "496", :c_ions "0",
:b_score "1", :y_ions "4", :a_score "0", :z_ions "0", :id "2693.1.1", :delta
"0.0024", :post "GLTG", :a_ions "0", :seq "GADGAPGKDGVR", :end 148,
:z_score "0", :x_score "0", :y_score "1", :c_score "0"}
clj-tandem.core> 
```

The functions 'protein-report' and 'peptide-report' return a
collection of strings representing csv separated lines describing
peptide and protein identifications.

You can run X! Tandem with the 'xtandem' function which takes paths to
a FASTA database, a spectra file and a map of X! Tandem
parameters. There are a set of default parameters (which can be seen
by using the 'default-inputs' function). By default the parameters are
for the TPP version of X! Tandem i.e k-score scoring, but this can be
changed using the parameters map. By default the function assumes
there is a program 'tandem' in the PATH variable and uses this as the
executable but this can be changed using the :tandem keyword and
specifying the path to the X! Tandem executable desired. The function
returns the path to the output file (which is the spectra file with
the extension replaced with '.tandem.xml').

## License

Copyright © 2016 Jason Mulvenna

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
