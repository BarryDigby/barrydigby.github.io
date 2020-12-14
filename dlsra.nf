#!/usr/bin/env nextflow

Channel
    .fromSRA('SRP079185')
    .view()
