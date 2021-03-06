mkdir -p plots.trash
mv plots/*.pdf plots.trash
mkdir -p tables.trash
mv tables/* tables.trash
mkdir -p plots.an/results
set -e // Exit when any command fails
rp run/higgsino/plot_n_minus_1.exe --sample search --year run2 --unblind --unblind_signalregion
rp run/higgsino/plot_search_unblind.exe --year run2 --unblind --unblind_signalregion --string_options "paper_style"
rp run/higgsino/plot_kappas.exe --sample search --unblind --year run2 --scen data --string_option do_zbi
scons&& run/higgsino/plot_kappas.exe --sample search --year run2 --scen mc
mv plots/* plots.an/results
mv tables/*.tex plots.an/results
