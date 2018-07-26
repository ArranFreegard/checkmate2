import os, sys
import subprocess
import signal
import ConfigParser
from math import sqrt
from copy import deepcopy
from advprint import AdvPrint
from info import Info
from resultcollector import ResultCollector

# note that for some strange reason, I cannot import "events" in the header as this results in a circular reference to Info.

class Process:
    name = '' # User defined name
    
    have_xsect = False # cross section provided by user?
    xsec = 0.0 # Cross section
    xsec_unit = 'pb' # Unit for cross section, pb or fb
    
    have_xerr = False # cross section error provided by user?
    xerr = 0.0 # Error
    xerr_unit = 'pb' # Unit for cross section error. can also be '%'

    have_kfac = False # K factor provided by user?
    kfac = 0.0  # k-factor (to be used when using internal Pythia or Madgraph processes)

    eventsList = list()

    have_reweighting = False
    reweighting_spectrum_dir = ''
    n_reweighting_targets = 0
    have_matrix_element_reweighting = False
    reweighting_me_proc_name = ''
    reweighting_pdf_name = ''
    reweighting_writing_mode = False
    reweighting_detector_simulation = False
    
    result_output_file = ""

    # Question: Why are there static members of the class? Surely the name of a process does not need to be a property of the class but the object.
    # And that is also why the exact same variables are created again at the object level in __init__...

    def __init__(self, name):
        self.name = name
        self.have_xsect = False
        self.xsec = 0.0 # Cross section
        self.xsec_unit = 'pb' # Unit for cross section, pb or fb
        self.have_xerr = False
        self.xerr = 0.0 # Error
        self.xerr_unit = 'pb' # Unit for cross section error. can also be '%'
        self.have_kfac = False
        self.kfac = 0.0  # k-factor (to be used when using internal Pythia or Madgraph processes)
                
        self.eventsList = list()
            
        self.have_reweighting = False
        self.reweighting_spectrum_dir = ''
        self.n_reweighting_targets = 0
        self.have_matrix_element_reweighting = False
        self.reweighting_me_proc_name = ''
        self.reweighting_pdf_name = ''
        self.reweighting_writing_mode = False
        self.reweighting_detector_simulation = False

        self.result_output_file = os.path.join(Info.paths["output_evaluation"], self.name+"_processResults.txt")
    
    
    def combine_processes(self, other):
        """ In case of an "add" run, the user might add new events to a given process """
        if self.name != other.name:
            AdvPrint.cerr_exit("'add' feature tried to combine incompatible processes (different names)")
            
        for other_ev in other.eventsList:
            self.eventsList.append(other_ev)
    

    def printInfo(self):
        AdvPrint.cout("\tProcess Name: "+self.name)
        if self.have_xsect:
            AdvPrint.cout("\tInput Cross section: "+str(self.xsec)+" "+str(self.xsec_unit))
        if self.have_xerr:
            AdvPrint.cout("\tInput cross section error: "+str(self.xerr)+" "+str(self.xerr_unit))
        elif self.have_kfac:
            AdvPrint.cout("\tInput KFactor: "+str(self.kfac))
            
        eventsList = [ef for ef in self.eventsList]
        if len(eventsList) != 0:   
            AdvPrint.cout("\tAssociated event files and/or Monte-Carlo generation runs:")
            for ef in eventsList:
                ef.printInfo()
                AdvPrint.cout("")
                
            AdvPrint.cout("")

    def checkInputConsistency(self):
        """ Checks if the amount of input information given is sufficient to run CheckMATE2 """        
    
        # Case 1: no process, no eventfile, no setting       
        if len(self.eventsList) < 1:
            self.printInfo()
            AdvPrint.cerr_exit("\t "+self.name+"::checkInputConsistency()::  \n \t"
            "There are no events to generate and/or analyse")
            
        # Case 3: Both k-factor and cross section are non-zero
        if self.xsec != 0.0 and self.kfac != 0.0:
            AdvPrint.cerr_exit("\t "+self.name+":checkInputConsistency()::  \n \t "
            "Either enter Kfactor or total cross section (which might be Kfactor * LO-cross section)!")
    
        for e in self.eventsList:
            from events import DelphesEvents
            if isinstance(e, DelphesEvents) and len(Info.used_experiments) > 1:
                AdvPrint.cerr_exit("\t "+self.name+"::checkInputConsistency()::  \n \t"
                                   "If you provide a .root file you cannot run analyses that need different detector settings."
                                   "Please only run analyses of the experiment the .root file was simulated with "
                                   "in Delphes!")
                
        for events in self.eventsList:
            events.check()
        
    def prepare(self):
        for events in self.eventsList:
            if events.processed:
                continue
            events.prepare()
    
    def prepareFritzInputFile(self, config, events):
        config, secname = events.prepareFritzInputFile(config)
        # if cross section is provided, override config
        if self.have_xsect:                
            config.set(secname, "xsect", str(self.xsec*Info.unit(self.xsec_unit)))
        if self.have_xerr:
            if self.xerr_unit == "%":
                config.set(secname, "xsectErrFactor", str(self.xerr/100.))
            else:
                config.set(secname, "xsectErr", str(self.xerr*Info.unit(self.xerr_unit)))
        else: # fritz wants this, lets give it to him
            config.set(secname, "xsectErr", str(0))
            
        if self.have_kfac:
            config.set(secname, "kfactor", str(self.kfac))
                        
        return config


    def prepareFritzReweighting(self, config, events):    
        name = "ReweightingHandler: "+events.identifier
        config.add_section(name)
        
        config.set(name, "outputDirectory", Info.paths["output_reweighting"])
        settings = os.path.join(Info.paths["output_reweighting"],"reweighting_"+events.identifier+".ini")
        config.set(name, "settings", settings)
        logFile = os.path.join(
                Info.paths["output_reweighting"],
                "reweighting_"+events.identifier+".log"
                )
        config.set(name, "logFile", logFile)

        from events import Pythia8Events, DelphesEvents   
        if isinstance(events, Pythia8Events):
            config.set(name, "pythiaHandler", events.identifier)
        elif isinstance(events, DelphesEvents):
            AdvPrint.cerr_exit("Reweighting does not support ROOT input files yet.")
        else:
            config.set(name, "eventFile", events.identifier)

        if not os.path.isdir(self.reweighting_spectrum_dir):
            AdvPrint.cerr_exit("reweighting_spectrum_dir does not exist: "+self.reweighting_spectrum_dir)
        slha_files = os.listdir(self.reweighting_spectrum_dir)

        cross_section_info = dict()
        if os.path.exists(os.path.join(self.reweighting_spectrum_dir,"xsec.dat")):
            with open(os.path.join(self.reweighting_spectrum_dir,"xsec.dat"),"r") as xsec_file:
                for line in xsec_file:
                    target_name, xsec, xsec_err = line.strip().split("\t")
                    cross_section_info[target_name] = (xsec, xsec_err)
            slha_files.remove("xsec.dat")

        if any([os.path.splitext(slha_file)[1]!=".slha" for slha_file in slha_files]):
            AdvPrint.cerr_exit("One of the given spectrum files is not an .slha file.")
        if not "base.slha" in slha_files:
            AdvPrint.cerr_exit("No base.slha file given.")
        slha_files.remove("base.slha")

        if len(slha_files) < 1:
            AdvPrint.cerr_exit("No reweighting target .slha files given.")
        self.n_reweighting_targets = len(slha_files)
        
        config.set(name, "ntargets", str(self.n_reweighting_targets))

        reweightingConfig = ConfigParser.RawConfigParser()
        globalSectionName = "Section: Global"
        reweightingConfig.add_section(globalSectionName)
        reweightingConfig.set(globalSectionName, "matrix_element_reweighting", "true")
        reweightingConfig.set(globalSectionName, "detector_simulation", self.reweighting_detector_simulation)
        reweightingConfig.set(globalSectionName, "writing_mode", self.reweighting_writing_mode)
        reweightingConfig.set(globalSectionName, "pdf_name", self.reweighting_pdf_name)
        if self.have_matrix_element_reweighting:
            reweightingConfig.set(globalSectionName, "procname", self.reweighting_me_proc_name)
        
        slhaFilesSectionsName = "Section: SLHAFiles"
        reweightingConfig.add_section(slhaFilesSectionsName)
        reweightingConfig.set(slhaFilesSectionsName, "base", os.path.abspath(os.path.join(self.reweighting_spectrum_dir, "base.slha")))
        for itarget,slha_file in enumerate(slha_files):
            full_path = os.path.abspath(os.path.join(self.reweighting_spectrum_dir,slha_file))
            reweightingConfig.set(slhaFilesSectionsName, "target{}".format(itarget+1), full_path)

        if len(cross_section_info) > 0:
            crossSectionsSectionsName = "Section: CrossSections"
            crossSectionErrorsSectionsName = "Section: CrossSectionErrors"

            reweightingConfig.add_section(crossSectionsSectionsName)
            reweightingConfig.add_section(crossSectionErrorsSectionsName)

            for itarget, (target_name,(xsec,xsec_err)) in enumerate(cross_section_info.iteritems()):
                target_name = "base" if itarget==0 else "target{}".format(itarget)
                reweightingConfig.set(crossSectionsSectionsName, target_name, xsec)
                reweightingConfig.set(crossSectionErrorsSectionsName, target_name, xsec_err)

        with open(settings, 'wb') as outfile: 
            reweightingConfig.write(outfile)

        return config

    
    def prepareFritzDelphes(self, config, events):
        for atype in Info.used_experiments:
            name = "DelphesHandler: "+atype
            config.add_section(name)
            settings = Info.files["delphes_global_config"][atype]
            config.set(name, "settings", settings)
            logFile = os.path.join(
                    Info.paths["output_delphes"],
                    "delphes_"+events.identifier+".log"
                    )
            config.set(name, "logFile", logFile)
            if Info.flags["write_delphes_events"]:
                outputFile = os.path.join(
                        Info.paths["output_delphes"],
                        events.identifier+"_"+atype+".root"
                        )
                config.set(name, "outputFile", outputFile)
            
            if self.have_reweighting:
                config.set(name, "reweightingHandler", events.identifier)
            else:
                from events import Pythia8Events    
                if isinstance(events, Pythia8Events):
                    config.set(name, "pythiaHandler", events.identifier)
                else:
                    config.set(name, "eventFile", events.identifier)


    def prepareFritzAnalysisHandler(self, config, events):
        if Info.flags["skipanalysis"]:
            return
        for atype in Info.used_experiments:
            name = "AnalysisHandler: "+atype
            config.add_section(name)
            analysis_type=atype
            config.set(name, "analysistype", analysis_type)
            config.set(name, "outputPrefix", events.identifier)
            config.set(name, "outputDirectory", Info.paths["output_analysis"])
            config.set(name, "logFile", Info.files["analysis_log"])
            from events import DelphesEvents
            if isinstance(events, DelphesEvents):
                config.set(name, "eventFile", events.identifier)
            else:
                config.set(name, "delphesHandler", atype)

    def prepareFritz(self):
        for events in self.eventsList:
            if events.processed:
                continue
            globalconfig = Info.config
            # add process dependent global settings
            name = "Global"
            if events.maxEvents > 0:
                if not globalconfig.has_section(name):
                    globalconfig.add_section(name)
                globalconfig.set(name, "nevents", events.maxEvents)
            elif globalconfig.has_section("Global"):
                globalconfig.remove_option(name, "nevents") # unfortunately this globalconfig construction is badly designed. Each event has its own maxevents but there can be only one globalconfig, so each one has to make sure that globalconfig is set correctly
            config = ConfigParser.RawConfigParser()
            
            self.prepareFritzInputFile(config, events)
            self.prepareFritzReweighting(config, events)
            self.prepareFritzDelphes(config, events)
            self.prepareFritzAnalysisHandler(config, events)
            path = os.path.join(Info.paths["output_fritz"],events.identifier+".ini")

            with open(path, 'wb') as file: 
                globalconfig.write(file)
            with open(path, 'ab') as file:
                config.write(file)
            events.configFile = path
    
    def runFritz(self):
        self.prepareFritz()
        """ Runs Fritz """
        from events import MG5Events
        for event in self.eventsList:
            if event.processed:
                continue
            fritz_command = Info.files["fritz_bin"]+" "+event.configFile
            result = subprocess.Popen(fritz_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #result = subprocess.Popen("echo foobar", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            maxlen = 0
            try:
                for line in iter(result.stdout.readline, b''):
                    # Print to logfile. If it does not start with the Fritz::Global prefix, it comes from MG5 and should be redirected
                    AdvPrint.mute()
                    if not line.startswith("|~| ") and isinstance(event, MG5Events):
                        AdvPrint.set_cout_file(os.path.join(Info.paths["output_mg5"], "mg5amcatnlo_"+event.identifier+".log"))
                    elif "PYTHIA Rndm::dumpState" in line:
                        AdvPrint.set_cout_file(os.path.join(Info.paths["output_pythia"], "pythia_"+event.identifier+".log"))
                    else:
                        line = line.replace("|~| ", "")
                        AdvPrint.set_cout_file(os.path.join(Info.paths["output_fritz"], "fritz_"+event.identifier+".log"))
                    AdvPrint.cout(line.rstrip())
                    AdvPrint.set_cout_file("#None")
                    AdvPrint.unmute()
                    
                    # We should not exceed the terminal terminal_width:
                    terminal_width = AdvPrint.get_terminal_width()
                    print_line = " |-> "+str(line.strip())
                    print_line.replace("\t", "    ")
                    while "\r" in print_line:
                        print_line = print_line[print_line.index("\r"):]
                        
                    len_of_print_line = len(print_line)
                    maxlen = max(maxlen, len(print_line))                
                    
                    # As we print line by line in the same terminal row, we have to add spaces if the curr line is shorter than a line before
                    fill_spaces = ""
                    if len(print_line) < maxlen:
                        fill_spaces = " "*(maxlen-len(print_line))
                        
                    # if line is too long, make it shorter by appropriate amoung
                    if len(print_line+fill_spaces) >= terminal_width and len(print_line) <= terminal_width:
                        fill_spaces = " "*(terminal_width - len(print_line)-1)
                    elif len(print_line) > terminal_width:
                        fill_spaces = ""
                        print_line = print_line[:terminal_width-4]+"..."                    
                        
                    AdvPrint.cout("\r"+print_line+fill_spaces+"\x1b[0m\r", "nlb")
            except KeyboardInterrupt:
                AdvPrint.cout("Caught Keyboard Signal. Aborting Fritz")
                result.send_signal(signal.SIGTERM)
            
            for line in iter(result.stderr.readline, b''):
                AdvPrint.unmute()
                AdvPrint.set_cout_file(Info.files['fritz_log'])
                # remove nasty ROOT6-CLING warnings from on-screen output 
                if "cling::AutoloadingVisitor::InsertIntoAutoloadingState:" in line:
                    AdvPrint.mute()
                elif "Missing FileEntry for ExRootAnalysis" in line:
                    AdvPrint.mute()
                elif "requested to autoload type" in line:
                    AdvPrint.mute()
                AdvPrint.cout(line.rstrip()+"")
                AdvPrint.set_cout_file("#None")
                AdvPrint.unmute()                
            AdvPrint.cout("")
            # Abort if there was an error
            result.wait()
            if result.returncode != 0:
                AdvPrint.cerr_exit("Fritz returned with error. Check logfiles in result folder for more information!")
            
            # Remove all empty analysisstdout files
            for f in [x for x in os.listdir(Info.paths['output_analysis']) if x.startswith("analysisstdout")]:
                if os.stat(os.path.join(Info.paths['output_analysis'], f)).st_size == 0:
                    os.remove(os.path.join(Info.paths['output_analysis'], f))

            # Associate result files to event
            for a in Info.analyses:
                for itarget in xrange(self.n_reweighting_targets+1):
                    a_key = a if itarget == 0 else "{}_target{}".format(a,itarget)

                    event.analysis_signal_files[a_key] = os.path.join(Info.paths['output_analysis'], event.identifier+'_'+a_key+'_signal.dat')
                    if os.path.isfile(event.analysis_signal_files[a_key]):
                        AdvPrint.format_columnated_file(event.analysis_signal_files[a_key])
                    
                    event.analysis_cutflow_files[a_key] = os.path.join(Info.paths['output_analysis'], event.identifier+'_'+a_key+'_cutflow.dat')
                    if os.path.isfile(event.analysis_cutflow_files[a_key]):
                        AdvPrint.format_columnated_file(event.analysis_cutflow_files[a_key])
            
            # finish
            event.processed = True

    def run(self):
        self.runFritz()
    
    def get_resultCollectors(self):        
        """ Gathers results from all events"""
        # setup resultCollector object
        resultCollectors_pr = dict()
                   
        # loop over all associated events and average results in all resultCollectors
        for event in self.eventsList:
            resultCollectors_ev = event.get_resultCollectors()

            for analysis, signal_regions in resultCollectors_ev.iteritems():                
                if not analysis in resultCollectors_pr:
                    resultCollectors_pr[analysis] = dict()

                for sr, rc in signal_regions.iteritems():
                    if not sr in resultCollectors_pr[analysis]:
                        resultCollectors_pr[analysis][sr] = ResultCollector(self.name, analysis, sr)

                    resultCollectors_pr[analysis][sr].add_and_average(rc)
                    
        
        # Write process file, if wanted
        if Info.parameters["ProcessResultFileColumns"] != []: 
            AdvPrint.mute()        
            AdvPrint.set_cout_file(self.result_output_file, True)
            for col in Info.parameters["ProcessResultFileColumns"]:
                AdvPrint.cout(col+"  ", "nlb")
            AdvPrint.cout("")
            for a in sorted(resultCollectors_pr.keys()):
                for sr in sorted(resultCollectors_pr[a].keys()):
                    AdvPrint.cout(resultCollectors_pr[a][sr].line_from_data(Info.parameters["ProcessResultFileColumns"]))
            AdvPrint.format_columnated_file(self.result_output_file)
            AdvPrint.set_cout_file("#None")
            AdvPrint.unmute()
                
        return resultCollectors_pr
    
