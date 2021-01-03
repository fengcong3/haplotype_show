function HaplotypeGraph(argsMap) {
    var self = this;

    var containerId;
    var container;
    let _init = function () {
        containerId = getRequiredVar(argsMap, "containerId");
        container = document.querySelector("#" + containerId);

        // margin[0] = getOptionalVar(argsMap, 'marginTop', 20) // marginTop allows fitting the actions, date and top of axis labels
        // margin[1] = getOptionalVar(argsMap, 'marginRight', 10)
        // margin[2] = getOptionalVar(argsMap, 'marginBottom', 10) // marginBottom allows fitting the legend along the bottom
        // margin[3] = getOptionalVar(argsMap, 'marginLeft', 30) // marginLeft allows fitting the axis labels
    }

    _init();
}