FROM 4geniac/almalinux:9.3_conda-py311_23.10.0-1

LABEL gitUrl="https://gitlab-ci-token:64_H4txhYSN_Gya4RJTC_Gr@gitlab.curie.fr/data-analysis/vegan.git"
LABEL gitCommit="9d62c26c025ae92eba1fef22186357f0203e93fb/devel"

ENV R_LIBS_USER "-"
ENV R_PROFILE_USER "-"
ENV R_ENVIRON_USER "-"
ENV PYTHONNOUSERSITE 1
ENV PATH $PATH
ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
ENV BASH_ENV /opt/etc/bashrc


# real path from projectDir: transition.yml
ADD transition.yml /opt/transition.yml

RUN dnf clean all \
&& CONDA_ROOT=$(conda info --system | grep CONDA_ROOT | awk '{print $2}') \
&& micromamba env create --root-prefix ${CONDA_ROOT} -f /opt/transition.yml \
&& echo -e "#! /bin/bash\n\n# script to activate the conda environment transition" > ~/.bashrc \
&& echo "export PS1='Docker> '" >> ~/.bashrc \
&& conda init bash \
&& echo "conda activate transition" >> ~/.bashrc \
&& mkdir -p /opt/etc \
&& cp ~/.bashrc /opt/etc/bashrc \
&& conda clean -y -a \
&& micromamba clean -y -a 

ENV PATH /usr/local/conda/envs/transition/bin:$PATH

