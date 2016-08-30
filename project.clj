(defproject clj-tandem "0.1.6"
  :description "Parser for X! Tandem XML formatted result files."
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[me.raynes/fs "1.4.6"]
                 [org.clojars.hozumi/clj-commons-exec "1.2.0"]
                 [org.clojure/clojure "1.8.0"]
                 [org.clojure/data.xml "0.0.8"]
                 [org.clojure/data.zip "0.1.1"]
                 [org.clojars.hozumi/clj-commons-exec "1.2.0"]]
  :jvm-opts ["-Xmx2000m"])
