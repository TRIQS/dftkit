# Command line scripts.  install(PROGRAMS) sets the execute bit, unlike the
# install(FILES) that dft_tools used for init_dmftpr.
install(PROGRAMS init_dmftpr DESTINATION bin)
