;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "tikzlibraryquantikz2.code"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("xargs" "") ("ifthen" "") ("xstring" "") ("xparse" "") ("etoolbox" "") ("mathtools" "") ("pgfmath" "") ("environ" "")))
   (TeX-run-style-hooks
    "xargs"
    "ifthen"
    "xstring"
    "xparse"
    "etoolbox"
    "mathtools"
    "pgfmath"
    "environ")
   (TeX-add-symbols
    "endtikzcd"
    '("gategroup" ["argument"] 1)
    '("wave" ["argument"] 0)
    '("gateoutput" ["argument"] 1)
    '("gateinput" ["argument"] 1)
    '("midstick" ["argument"] 1)
    '("rstick" ["argument"] 1)
    '("lstick" ["argument"] 1)
    '("slice" ["argument"] 1)
    '("qwbundle" ["argument"] 1)
    '("setwiretype" ["argument"] 1)
    '("ghost" ["Text"] ["Text"] "Text")
    '("gate" ["Text"] ["Text"] ["Text"] "Text")
    '("makeebit" ["Text"] "Text")
    '("swap" ["Text"] "Text")
    '("expandedwire" "Text" "Text" "Text" "Text")
    '("octrl" ["Text"] "Text")
    '("ctrl" ["Text"] "Text")
    '("ground" ["Text"])
    '("trash" ["Text"])
    '("measure" ["Text"])
    '("inputD" ["Text"])
    '("meterD" ["Text"])
    '("measuretab" ["Text"])
    '("metercw" ["Text"] ["Text"] ["Text"] "Text")
    '("meter" ["Text"] ["Text"] ["Text"] "Text")
    '("targX" ["Text"] "Text")
    '("targ" ["Text"] "Text")
    '("ocontrol" ["Text"] "Text")
    '("control" ["Text"] "Text")
    '("phase" ["Text"] "Text")
    '("cwbend" "Text")
    '("wireoverride" ["Text"] ["Text"] "Text")
    '("wire" ["Text"] ["Text"] ["Text"] "Text")
    '("myrepeat" ["Text"] "Text" "Text")
    '("vqwexplicit" 2)
    '("braket" 2)
    '("proj" 1)
    '("bra" 1)
    '("ket" 1)
    '("setmiddle" 1)
    '("ifcsstringeitheror" 4)
    '("permute" 1)
    '("push" 1)
    '("hphantomgate" 1)
    '("phantomgate" 1)
    '("pgfdeclareanchoralias" 3)
    '("vqw" 1)
    '("vcw" 1)
    '("importwiretypes" 1)
    '("fullwire" 4)
    "resetwiretypes"
    "qw"
    "cw"
    "linethrough"
    "wirepair"
    "gate"
    "slice"
    "make"
    "sliceallr"
    "sliceallvr"
    "groupinput"
    "groupoutput"
    "mginput"
    "mgoutput"
    "wave"
    "gategroup"
    "DivideRowsCols"
    "forceredefine"
    "tikzcdmatrixname"
    "errmessage"
    "pgfaddtoshape"
    "anchor"
    "firstpoint"
    "secondpoint"
    "verticaltext"
    "xvvv"
    "tmp"
    "xxvvv"
    "ifnodedefined"
    "firstlist"
    "quantwires"
    "background"
    "toswap"
    "DisableMinSize"
    "leftbrace"
    "rightbrace"
    "resetstyles"
    "maketransparent")
   (LaTeX-add-environments
    "tikzpicture"
    '("quantikz" LaTeX-env-args ["Text"]))
   (LaTeX-add-counters
    "dummy"
    "aaa"
    "wirenumberpermute")
   (LaTeX-add-lengths
    "myl"
    "myh"
    "myd"))
 :latex)

