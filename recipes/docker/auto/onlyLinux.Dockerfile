FROM 4geniac/almalinux:9.3

LABEL gitUrl="https://gitlab-ci-token:64_H4txhYSN_Gya4RJTC_Gr@gitlab.curie.fr/data-analysis/vegan.git"
LABEL gitCommit="9d62c26c025ae92eba1fef22186357f0203e93fb/devel"

ENV R_LIBS_USER "-"
ENV R_PROFILE_USER "-"
ENV R_ENVIRON_USER "-"
ENV PYTHONNOUSERSITE 1
ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
