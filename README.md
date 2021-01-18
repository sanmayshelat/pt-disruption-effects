# README #

Modelling public transport disruption propagation effects using: 

* A fast, mesoscopic assignment model that accounts for on-board crowding and denied boarding.
* A new criticality indicator called spatial criticality

See abstract below for more details about the study. If you use this, please cite: Shelat, S., Cats, O., 2017. Measuring spill-over effects of disruptions in public transport networks, Models and Technologies for Intelligent Transportation Systems (MT-ITS), 5th IEEE International Conference on. IEEE, pp. 756-761.

### Data Required ###

Expects transit (OD - station level) demand and supply data. Sample datafiles available upon request to s.shelat@tudelft.nl


### How to: ###

* Create datafiles as required.
* Run the main file.


### Abstract ###

Transit is vital to the functioning of most major cities in the world today and therefore evaluation of the robustness - capacity to withstand disruptions with minimal impact to the system - of transit networks is essential. In this paper, we study the spatial extent of link disruption impacts in urban public transport networks (PTNs) due to spill-over effects that occur as affected passengers choose alternative routes. Quantifying link criticality in terms of its spatial propagation effects is important for prioritizing robustness measures as well as devising methods to encapsulate disruption spill-over effects. To this end, a new local criticality indicator - spatial criticality - that is based on: (i) the magnitude of relative change in load on other links due to closure of a link and (ii) the topological distance between these other links and the disrupted link is introduced. Further, a stochastic user equilibrium, PTN assignment model is developed. The model explicitly considers the service layer of the PTNs and accounts for on-board crowding and denied boarding at stops to represent passenger spill-over effects. Finally, to demonstrate the proposed indicator and investigate its relationship with conventional, local topological indicators, the urban rail-bound network of Amsterdam is used as a case study.