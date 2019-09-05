# pfc
## Tools for phase field crystal modeling of two-dimensional materials.

This repository contains various tools for phase field crystal modeling of two-dimensional materials; see my thesis (http://urn.fi/URN:ISBN:978-952-60-8608-8), our manuscript on heterostructures (https://arxiv.org/abs/1908.05564) and the references therein. These tools are cleaned-up versions of the many, many tools I have written and used throughout my doctoral work on the topic. Some examples include: pfc-init.jar for constructing initial states, pfc-relax for relaxing the systems and coordinator.jar for converting relaxed PFC density fields into atomic coordinates. Various visualization tools are also included. To use all of the tools, Java, MPI, FFTW3 and POV-Ray are required. Each tool prints out instructions if run without input arguments. Two demonstrations of how to use the tools are also provided. One considers graphene and the other lateral heterostructures of graphene and hexagonal boron nitride. Simply run the scripts provided (honey.sh and hetero.sh).

Note that I am now leaving academia and had limited time to prepare these tools for publication. Expect bugs.
