;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "paper_body_v1"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("revtex4-2" "aps" "prx" "showpacs" "preprintnumbers" "twocolumn" "superscriptaddress" "10pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("packages" "")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "setfloatlink")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "homepage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "email")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "../visual_elements/figs/bipartite_graph"
    "../visual_elements/figs/stabilizer_diagrams"
    "revtex4-2"
    "revtex4-210"
    "packages")
   (TeX-add-symbols
    '("includesupp" 1)
    '("includemain" 1)
    "id"
    "Stot"
    "Stota"
    "StotaN"
    "pc"
    "ps"
    "mpipks"
    "julich"
    "regensburg"
    "thinline")
   (LaTeX-add-labels
    "fig:model"
    "sec:model"
    "eq:Stot"
    "eq:normalization"
    "fig:combined"
    "sec:random"
    "eq:random_policy"
    "sec:greedy"
    "eq:greedy_policy"
    "sec:learned_str"
    "sec:rl_circuit_control"
    "fig:rl_scheme"
    "sec:informed_strategies"
    "fig:rl_fullinfo"
    "eq:prob_relation"
    "sec:partial_info"
    "fig:rl_p"
    "sec:discussion_outlook"
    "app:captions"
    "fig:example_snapshot"
    "app:nonequilibrium"
    "app:supp_results"
    "fig:transient"
    "app:enhanced_transitions"
    "fig:extended_results"
    "app:binder"
    "fig:greedy_binder"
    "app:pyramid_strategies"
    "fig:pyramids_human"
    "app:vanishing_info"
    "fig:vanishing_info"
    "app:rl"
    "app:rl_environment"
    "alg:step"
    "app:rl_ppo"
    "app:rl_policy"
    "fig:encoder"
    "app:rl_training"
    "app:comp_resources"
    "app:rl_hyperparameters"
    "tab:hyperparams"
    "app:stabilizer_circuits"
    "eq:Pauli_generating_set"
    "th:aaronson"
    "app:simulating_stabilizers"
    "tab:operators_encoding"
    "eq:identity-tableau"
    "tab:cnot_example"
    "app:disentangling_cliffords"
    "eq:clipped_condition"
    "app:clipping_algorithm")
   (LaTeX-add-environments
    "theorem")
   (LaTeX-add-bibliographies
    "references"))
 :latex)

