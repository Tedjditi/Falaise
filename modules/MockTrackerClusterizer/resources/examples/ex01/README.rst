MockTrackerClusterizer Falaise module
=====================================

Introduction
------------

This  example shows  how to  build a  ``flreconstruct`` pipeline  that
processes the  SuperNEMO simulated data  (SD bank from  the flsimulate
output)  to produce  mock calibration  data  (CD bank)  then the  mock
tracker clusterizer (TCD bank).

Configuration file(s):

  * ``config/pipeline.conf``   :    the   ``flreconstruct``   pipeline
    configuration  files that  describes the  chain of  modules to  be
    processed.

Usage
-----

  1. First setup the Falaise software.

  2. Check the geometry of the SuperNEMO demonstrator (not mandatory): ::

      $ bxgeomtools_inspector \
          --datatools::logging "warning" \
          --load-dll Falaise@${FALAISE_INSTALL_DIR}/lib/x86_64-linux-gnu \
          --datatools::resource-path "falaise@${FALAISE_INSTALL_DIR}/share/Falaise-1.0.0/resources" \
          --manager-config "@falaise:config/snemo/demonstrator/geometry/3.0/manager.conf"

  3. Generate a set of simulated events: ::

      $ flsimulate \
        --experiment "demonstrator" \
        --vertex-generator "source_strips_bulk" \
        --event-generator  "Se82.0nubb" \
	--number 10 \
	--output-file "Se82_0nubb-source_strips_bulk_SD.brio"

  4. Run the reconstruction pipeline:

      * Browse the manual for the mock tracker clustering module: ::

          $ bxocd_manual \
            --load-dll Falaise@${FALAISE_INSTALL_DIR}/lib/x86_64-linux-gnu \
            --load-dll Falaise_MockTrackerClusterizer@${FALAISE_INSTALL_DIR}/lib/x86_64-linux-gnu/Falaise/modules \
            --action show snemo::reconstruction::mock_tracker_clustering_module > falaise_mtc.rst
          $ pandoc -w html falaise_mtc.rst -o falaise_mtc.html
	  $ xdg-open falaise_mtc.html &


      * Run the reconstruction pipeline: ::

          $ flreconstruct \
            --input-file "Se82_0nubb-source_strips_bulk_SD.brio" \
	    --pipeline "config/pipeline.conf" \
            --output-file "Se82_0nubb-source_strips_bulk_SD-CD-TCD.brio"

        The  ``flreconstruct``  executable  should  perform  the  mock
        pre-reconstruction and then print  the structure of the processed
        event records  in the  terminal. You  should be  able to  see the
        ``"SD"``, ``"CD"`` and ``"TCD"`` banks.

        The output  file ``Se82_0nubb-source_strips_bulk_SD-CD-TCD.brio``
        is  thus ready  for  further reconstruction  (track fitting)  and
        analysis.
