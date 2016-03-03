(ns clj-tandem.core
  (:require [clojure.data.xml :refer [parse]]
            [clojure.data.zip.xml :refer [xml-> xml1-> text attr= attr]]
            [clojure.zip :refer [xml-zip node]]
            [clojure.string :refer [split]]))

(defn group-seq
  "Takes a buffered reader on a X! Tandem XML formatted file and
  returns a lazy list of zippers corresponding to models in the file."
  [reader]
  (->> (filter #(= (:tag %) :group) (:content (parse reader)))
       (map xml-zip)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; peptides and accessors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- protein-info
  "Returns a hash map of protein info from a peptide
  zipper. Including :id, :expect, :description etc."
  [pz]
  (let [[acc desc] (split (xml1-> pz :note (attr= :label "description") text)
                          #"\s+"
                          2)]
    (merge (->> (xml1-> pz node) :attrs)
           {:description desc :accession acc})))

(defn- peptide-info
  "Returns a hash map of peptide
  info. Including :id, :start, :end, :seq etc. Includes a :mod key
  which is a list of maps containing modification
  information (:type, :at and :modified) and a :peptide-info key
  containing protein level information."
  [pz]
  (let [pi (protein-info pz)
        pepi (->> (xml1-> pz :peptide :domain)
                  node
                  :attrs)]
    (merge pepi
           {:mods (->> (xml-> pz :peptide :domain :aa node)
                       (map :attrs))
            :rank (-> (:id pi)
                      (split #"\.")
                      second
                      Integer/parseInt)
            :start (Integer/parseInt (:start pepi))
            :end (Integer/parseInt (:end pepi))
            :expect (Float/parseFloat (:expect pepi))
            :protein-length (->> (xml1-> pz :peptide (attr :end))
                                 Integer/parseInt)
            :protein-info pi})))

(defn peptide-seq
  "Takes a group zipper and returns a collection of maps describing
  the peptides in a model."
  [gz]
  (->> (xml-> gz :protein)
       (map peptide-info)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; peptide report
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- peptide-headers
  []
  (->> '("Experiment_name" "Protein_accession" "Protein_description"
         "Peptide_ID" "Peptide_expect" "Precursor_mh"
         "Delta" "Rank" "Peptide_sequence" "Pre_sequence"
         "Post_sequence" "Peptide_start_index" "Peptide_end_index"
         "Peptide_mods" "Also_attributed")
       (interpose ",")
       (apply str)))

(def pep-rep-keys [:en :accession :description :id :expect :mh :delta :rank :seq :pre :post :start :end :mods :alsos])

(defn- peptide-collate
  [gz en]
  (let [accessions (->> (peptide-seq gz) (map #(:accession (:protein-info %))))]
    (->> (peptide-seq gz)
         (map #(merge (select-keys % pep-rep-keys)
                      {:en en
                       :accession (:accession (:protein-info %))
                       :description (str "\"" (:description (:protein-info %)) "\"")
                       :mods (if-let [m (->> (:mods %)
                                             (map (fn [x]
                                                    (str (:modified x)
                                                         "@"
                                                         (:type x)
                                                         (:at x))))
                                             seq)]
                               (apply str (interpose ";" m))
                               " ")
                       :protein-length (:protein-length %)
                       :protein-expect (:expect (:protein-info %))
                       :alsos (if-let [a (->> (remove (fn [x]
                                                        (= (:accession
                                                            (:protein-info %)) x))
                                                      accessions)
                                              seq)]
                                (->> (interpose ";" a) (apply str))
                                " ")})))))

(defn- get-peptides
  [reader en]
  (->> (group-seq reader)
       (mapcat #(peptide-collate % en))
       (group-by :accession)
       vals))

(defn peptide-report
  "Returns a collection of strings describing peptides. Only reports
  peptides where the identified protein has at least one peptide
  ranked 1."
  ([reader] (peptide-report reader "Experiment"))
  ([reader exp-name]
   (cons (peptide-headers)
         (->> (get-peptides reader exp-name)
              (filter #(some (fn [x] (= 1 (:rank x))) %))
              (sort-by count >)
              flatten
              (map #(->> (map (fn [x] (x %)) pep-rep-keys)
                         (interpose ",")
                         (apply str)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; protein report
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn protein-headers
  []
  (->> '("Experiment name" "Accession" "Description"
         "Expect" "Unique peptides" "Total peptides" "Cover")
       (interpose ",")
       (apply str)))

(defn- cover
  [up pl]
  (float (* 100
            (/ (->> (mapcat #(range (:start %) (+ 1 (:end %))) up)
                    set
                    count)
               pl))))

(def prot-rep-keys [:en :accession :description :expect :up :tp :cover])

(defn- unique-peps
  [peps]
  (let [t (group-by :seq (seq peps))]
    (->> (group-by :seq (seq peps))
         vals
         (map #(sort-by :expect < %))
         first)))

(defn- collate-proteins
  [peps en]
  (let [pep-count (count peps)
        up (unique-peps peps)
        c (cover up (:protein-length (first up)))]
    {:expect (:protein-expect (first up))
     :en en
     :accession (:accession (first up))
     :description (:description (first up))
     :up (count up)
     :tp pep-count
     :cover c}))

(defn protein-report
  "Returns a collection of strings describing protein identifications.
  Only reports peptides where the identified protein has at least one
  peptide ranked 1."
  ([reader] (protein-report reader "Experiment"))
  ([reader exp-name]
   (cons (protein-headers)
         (->> (get-peptides reader exp-name)
              (filter #(some (fn [x] (= 1 (:rank x))) %))
              (map #(collate-proteins % exp-name))
              (sort-by :up >)
              (map #(->> (map (fn [x] (x %)) prot-rep-keys)
                         (interpose ",")
                         (apply str)))))))
