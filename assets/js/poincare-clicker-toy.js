(function(root) {

  "use strict";

  function bMax(energy) {
    return Math.acos(-2.-energy);
  };

  function lbMax(energy) {
    return Math.sqrt(6.+2.*energy);
  };

  function lbSeparatrix(energy, b) {
    return Math.sqrt(2.*(2. + energy + Math.cos(b)));
  };

  function discrim(energy, b, lb) {
    return (4. + 2.*energy - lb*lb + 2.*Math.cos(b))*(3. - Math.cos(2.*b)) / 2.;
  };

  function insideSep(energy, b, lb) {
    return (discrim(energy, b, lb) >= 0);
  };

  function initConds(energy, b, lb) {
    // make sure we can get physical init conds
    var d = discrim(energy, b, lb);
    if (!(d >= 0.))
      throw 'initConds: b and lb are out of bounds for this energy';
    return [0., // a
            lb + lb*Math.cos(b) + Math.sqrt( d ), // la
            b, lb];
  };

  function ODERHS(t, xs) {
    var a = xs[0], la = xs[1], b = xs[2], lb = xs[3];

    var denom = 3. - Math.cos(2.* b);
    var cosb = Math.cos(b);
    var onePlusCosb = 1. + cosb;
    var threePlusTwoCosb = 3. + 2.* cosb;

    var adot = 2.*(la - lb*onePlusCosb) / denom;
    var bdot = 2.*(-la*onePlusCosb + lb*threePlusTwoCosb)/denom;
    var ladot = -2. * Math.sin(a) - Math.sin(a+b);
    var lbdot = -Math.sin(a+b) - 2. * lb * Math.sin(b) * (la - lb) / denom
        + 2. * (la*la - 2.*la*lb*onePlusCosb + lb*lb*threePlusTwoCosb) * Math.sin(2.*b) / denom / denom;

    return [adot, ladot, bdot, lbdot];
  };

  // Returns the x value where y crosses zero, using linear interpolation.
  function zeroCrossing(x0, y0, x1, y1) {
    var m = (y1-y0)/(x1-x0);
    var b = y1-m*x1;
    return -b/m;
  }

  /* Assume that y(x) is a linear function with values
   * y0 = y(x0), y1 = y(x1)
   * Then interpolate y(xprime)
  */
  function linearInterp(xprime, x0, y0, x1, y1) {
    var m = (y1-y0)/(x1-x0);
    var b = y1-m*x1;
    return m*xprime + b;
  }

  /* This is the heart of the toy.
   * Tries to yield npt _new_ points in the (b,lb) plane.
   * Returns an Array, hopefully of length npt+1, where the 0th
   * element is the initial point that was passed, and the rest are new.
   * Stops if maxSteps is reached, so it may yield fewer points.
   */
  function reapPoincarePoints(energy, npt, initB, initLB, deltaT, maxSteps) {
    var poincPoints = new Array();
    var foundPoints = 0, iter = 0;
    var curPoint = initConds(energy, initB, initLB), newPoint;
    var tCross, newB, newLB;

    // We will return the initial point as element 0
    poincPoints.push([initB, initLB]);

    // Fence-post issue with JXG.Math.Numerics.rungeKutta .
    // To get it to take one step, we need to tell it we want a total
    // of two steps on the interval, because it counts the 0th step.
    deltaT = 2.*deltaT;

    do {
      // Take an RK step
      newPoint = JXG.Math.Numerics.rungeKutta(
        'rk4',        // Butcher table
        curPoint,     // Initial conditions
        [0., deltaT], // time interval
        2,            // how many points
        ODERHS        // the RHS of the system
      );

      // Just get the new point
      newPoint = newPoint[1];

      // Check if there has been a zero-crossing for a, in the
      // positive direction
      if ((curPoint[0] < 0) && (newPoint[0] >= 0) ) {
        // We found a new point on the section!
        foundPoints++;

        // Find the approximate time of crossing
        tCross = zeroCrossing(0., curPoint[0], deltaT, newPoint[0]);

        // Interpolate the values of b, lb
        newB  = linearInterp(tCross, 0., curPoint[2], deltaT, newPoint[2]);
        newLB = linearInterp(tCross, 0., curPoint[3], deltaT, newPoint[3]);

        // store it
        poincPoints.push([newB, newLB]);
      };

      // Next
      iter++;
      curPoint = newPoint;

    } while ( (iter <= maxSteps) && (foundPoints < npt) );

    return poincPoints;
  };

  ////////////////////////////////////////////////////////////
  // UI related

  /* Emits the slider-drag-handler */
  function makeESliderDrag(controller) {
    return function() {
      controller.setenergy(controller.eslider.Value());
    };
  };

  /* Emits the slider-drag-handler */
  function makeNPtSliderDrag(controller) {
    return function() {
      controller.setnpt(controller.nptslider.Value());
    };
  };

  /* This emits a click handler for the Poincare section.
   * This is the UI-specific code.  It dispatches to
   * controller.handleTouch after translating to logical (b, lb) coordinates.
   * handleTouch needs to know if there is already a point there, and
   * if so, which one
   */
  function makePoincTouch(controller) {
    // See https://jsxgraph.uni-bayreuth.de/wiki/index.php/Browser_event_and_coordinates

    var board = controller.poincbox;
    var getMouseCoords = function(e, i) {
      var cPos = board.getCoordsTopLeftCorner(e, i),
          absPos = JXG.getPosition(e, i),
          dx = absPos[0]-cPos[0],
          dy = absPos[1]-cPos[1];

      return new JXG.Coords(JXG.COORDS_BY_SCREEN, [dx, dy], board);
    };

    return function(e) {
      var ptExists = false, i, coords, thePoint;

      if (e[JXG.touchProperty]) {
        // Is this right?
        if (e.touches.length > 1) // multi-touch, pass through to treat as pan or zoom?
          return;

        // Calling preventDefault() so that we don't also generate a
        // click event, see
        // https://developer.mozilla.org/en-US/docs/Web/API/Touch_events/Supporting_both_TouchEvent_and_MouseEvent
        e.preventDefault();

        // index of the finger that is used to extract the coordinates
        i = 0;
      }
      coords = getMouseCoords(e, i);
      var b = coords.usrCoords[1];
      var lb = coords.usrCoords[2];

      for (var el in board.objects) {
        if(JXG.isPoint(board.objects[el]) && board.objects[el].hasPoint(coords.scrCoords[1], coords.scrCoords[2])) {
          ptExists = true;
          thePoint = board.objects[el];
          break;
        }
      }

      controller.handleTouch(b, lb, ptExists, thePoint);

    };
  };

  // Add the separatrix curves to the plot
  function createSeparatrixCurves(controller) {
    /* q will run from 0 to 1, map it to the correct range of b */
    var bOfq = function(q){ return bMax(controller.energy)*(-1. + 2. * q); };
    var sepA = controller.poincbox.create('curve',
                 [bOfq,
                  function(q){ return lbSeparatrix(controller.energy, bOfq(q));},
                  0., 1.]);
    var sepB = controller.poincbox.create('curve',
                 [bOfq,
                  function(q){ return -lbSeparatrix(controller.energy, bOfq(q));},
                  0., 1.]);

    // Make them unclickable
    sepA.hasPoint = function(){return false; };
    sepB.hasPoint = function(){return false; };
  };

  function makeZoomHandler(controller) {
    return function(){
      var box = controller.poincbox.stopSelectionMode();
      // bbox has the coordinates of the selection rectangle.
      // Attention: box[i].usrCoords have the form [1, x, y], i.e.
      // are homogeneous coordinates.
      // Make sure to order min, max x,y so that orientation is preserved
      var xMin = Math.min(box[0].usrCoords[1],box[1].usrCoords[1]),
          xMax = Math.max(box[0].usrCoords[1],box[1].usrCoords[1]),
          yMin = Math.min(box[0].usrCoords[2],box[1].usrCoords[2]),
          yMax = Math.max(box[0].usrCoords[2],box[1].usrCoords[2]);
      // Set a new bounding box
      controller.poincbox.setBoundingBox([xMin, yMax, xMax, yMin], false);
      controller.isZoom100 = false;
      controller.updateButtonAbility();
    };
  };

  function makeClearClickHandler(controller) {
    return function(){
      // $('body').addClass('waiting');
      // setTimeout(function() {
        controller.clearPoints();
      //   $('body').removeClass('waiting');
      // }, 1);
    };
  };

  function makeMoreClickHandler(controller) {
    return function(){
      controller.morePointsFromLast();
    };
  };

  function makeUndoClickHandler(controller) {
    return function(){
      controller.popLastOrbit();
    };
  };

  function makeRedoClickHandler(controller) {
    return function(){
      controller.restoreOrbit();
    };
  };

  function makeZoom100ClickHandler(controller) {
    return function(){
      controller.initialZoom();
    };
  };

  function makePointOverHandler(controller, groupId, groupCSSClass) {
    var selectorText = "." + groupCSSClass;
    return function(e){
      var rules = controller.styleSheet.cssRules;
      for(var i=0; i<rules.length; i++) {
        // Find the correct rule in the stylesheet
        if(rules[i].selectorText == selectorText) {
          rules[i].style['fill'] = controller.extraStyles.hiliteColor;
          rules[i].style['rx'] = controller.extraStyles.hiliteSize;
          rules[i].style['ry'] = controller.extraStyles.hiliteSize;
        };
      };
    };
  };

  function makePointOutHandler(controller, groupId, groupCSSClass) {
    var selectorText = "." + groupCSSClass;
    return function(e){
      var rules = controller.styleSheet.cssRules;
      for(var i=0; i<rules.length; i++) {
        // Find the correct rule in the stylesheet
        if(rules[i].selectorText == selectorText) {
          rules[i].style['fill'] = controller.colorPal[controller.groupColMap[groupId]];
          rules[i].style['rx'] = controller.extraStyles.normalSize;
          rules[i].style['ry'] = controller.extraStyles.normalSize;
        };
      };
    };
  };

  function makeKeyHandler(controller) {
    return function(event) {
      return controller.handleKey(event);
    };
  }

  ////////////////////////////////////////////////////////////
  // Controller class

  // Constructor
  function PoincareClickerController(ctrlsboxName,buttonboxName,poincboxName) {

    if (!(this instanceof PoincareClickerController)) {
      return new PoincareClickerController(ctrlsboxName,buttonboxName,poincboxName);
    }

    this.setupBoxes(ctrlsboxName,buttonboxName,poincboxName);

  };

  // Controller prototype
  PoincareClickerController.prototype = {
    /* UI objects for the controls */
    ctrlsbox: {},
    eslider: {},
    nptslider: {},
    buttonbox: {},
    clearButton: {},
    moreButton: {},
    undoButton: {},
    redoButton: {},
    zoomButton: {},
    keyHelpButton: {},

    /* UI objects for the Poincare section box */
    poincbox: {},
    basePoincOpts: {
      boundingbox: [-2, 2, 2, -2],
      keepaspectratio: true,
      axis: false,
      grid: true,
      // renderer: 'canvas', // SVG seems to work better than canvas
      // pan: {enabled: true},
      showNavigation: true,
      showCopyright:  false},
    stylesheet: {},

    /* Control variables for making Poincare sections */
    energy: -1.9,
    npt: 250,
    deltaT: 0.03,
    maxSteps: 100000,

    /* Parameters for points during sliding energy */
    nSlidePts: 100,
    slideb: 0.01,
    slidelb: 0.01,

    isZoom100: true,

    /* styles */
    basePointStyle: {size: 0.5, sizeUnit: 'screen',
                       strokeWidth: 0,
                       color: '#000000',
                       fixed: true,
                       showInfobox: false,
                       name: '', withLabel: false},
    extraStyles: {normalColor: "#000000", normalSize: "1px",
                    hiliteColor: "#00bb00", hiliteSize: "1.5px",
                    defaultRuleString: ""},
    colorPal: ['#000000', '#4477aa', '#66ccee', '#228833', '#ccbb44', '#ee6677', '#aa3377'],

    /* Storage of points on Poincare section */
    groupCounter: 0,
    pointGroupList: new Array(),
    undonePointGroupList: new Array(),

    groupColMap: {},
    rotateColor: {},

    /* Public member functions */
    setupBoxes: {},
    setupCtrls: {},
    setupPoinc: {},

    setenergy: {},
    setnpt: {},

    handleTouch: {},
    handleKey: {},

    clearPoints: {},
    clearPointThresh: 500,
    _clearPointsSanely: {},
    _clearPointsViaFreeBoard: {},

    morePointsFromLast: {},
    popLastOrbit: {},
    restoreOrbit: {},

    bbox: {},
    padFactor: 1.05,
    initialZoom: {},

    updateButtonAbility: {},
  };

  // TODO Maybe setupBoxes should not be public
  PoincareClickerController.prototype.setupBoxes = function(ctrlsboxName,buttonboxName,poincboxName) {

    this.setupCtrls(ctrlsboxName,buttonboxName);

    this.setupPoinc(poincboxName);

    // Add keyboard handler
    document.body.addEventListener("keydown", makeKeyHandler(this));
  };

  PoincareClickerController.prototype.setupCtrls = function(ctrlsboxName,buttonboxName) {
    this.ctrlsbox = JXG.JSXGraph.initBoard(ctrlsboxName,
                                 {boundingbox:[0.,1.,1.,0.],
                                  axis:false,
                                  pan: {enabled: false},
                                  showNavigation: false,
                                  showCopyright:  false});
    this.ctrlsbox.suspendUpdate();

    this.eslider = this.ctrlsbox.create(
      'slider',
      [[0.05,.66],[0.7,.66],
       [-2.999,-1.909,-1.001]],
      {name: 'E', precision:3});

    this.nptslider = this.ctrlsbox.create(
      'slider',
      [[0.05,.33],[0.7,.33],
       [100,250,500]],
      {name: '# of points', snapWidth:1, precision:0});

    this.eslider.on('drag', makeESliderDrag(this));
    this.nptslider.on('drag', makeNPtSliderDrag(this));

    this.ctrlsbox.unsuspendUpdate();

    //////////////////////////////
    // Buttons

    this.buttonbox = document.getElementById(buttonboxName);

    this.clearButton = this.buttonbox.querySelector('#clear');
    this.moreButton = this.buttonbox.querySelector('#more');
    this.undoButton = this.buttonbox.querySelector('#undo');
    this.redoButton = this.buttonbox.querySelector('#redo');
    this.zoomButton = this.buttonbox.querySelector('#zoom100');
    this.keyHelpButton = this.buttonbox.querySelector('#keyHelp');

    this.clearButton
      .addEventListener('click',
                        makeClearClickHandler(this));

    this.moreButton
      .addEventListener('click',
                        makeMoreClickHandler(this));

    this.undoButton
      .addEventListener('click',
                        makeUndoClickHandler(this));

    this.redoButton
      .addEventListener('click',
                        makeRedoClickHandler(this));

    this.zoomButton
      .addEventListener('click',
                        makeZoom100ClickHandler(this));

    this.keyHelpButton
      .addEventListener('click',
                        keyHelp);
  };

  PoincareClickerController.prototype.setupPoinc = function(poincboxName) {

    this.poincbox = JXG.JSXGraph.initBoard(poincboxName, this.basePoincOpts);

    this.poincbox.suspendUpdate();
    this.poincbox.setBoundingBox(this.bbox(), false);

    var baxis = this.poincbox.create('axis', [[-4, 0], [4,0]],
        { name:'β',
          withLabel: true,
          ticks: {minorTicks:1, majorHeight:10, minorHeight:4},
          label: { position: 'rt',
                   offset: [-25, 20], }
        });
    var lbaxis = this.poincbox.create('axis', [[0, -4], [0,4]],
        { name:'ℓ_β',
          withLabel: true,
          ticks: {minorTicks:1, majorHeight:10, minorHeight:4},
          label: { position: 'rt',
                   offset: [-20, -25], }
        });

    // Is this the right way?  Make them non-clickable
    baxis.hasPoint = function(){return false; };
    lbaxis.hasPoint = function(){return false; };

    // For some reason JSXGraph creates a point at the origin, which
    // will capture our clicks there... this finds the point and
    // unsets its knowledge of owning (0,0)
    for (var el in this.poincbox.objects) {
      if(JXG.isPoint(this.poincbox.objects[el])
         && this.poincbox.objects[el].coords.usrCoords[1] == 0
         && this.poincbox.objects[el].coords.usrCoords[2] == 0)
      {
        this.poincbox.objects[el].hasPoint = function(){return false; };
      };
    };

    // For selection-zooming, see example at
    // https://jsxgraph.org/docs/symbols/JXG.Board.html#startSelectionMode
    this.poincbox.on('stopselecting', makeZoomHandler(this));

    createSeparatrixCurves(this);

    this.poincbox.on('down', makePoincTouch(this));

    // Add <style> element to the div containing the Poincare section
    // That is where we'll put the styles that control our point
    // groups for easy highlighting
    // See https://developer.mozilla.org/en-US/docs/Web/API/CSSStyleSheet/insertRule
    var styleEl = document.createElement('style');
    var poincDiv = document.getElementById(poincboxName);
    poincDiv.appendChild(styleEl);
    this.styleSheet = styleEl.sheet;
    this.basePointStyle.color = this.extraStyles.normalColor;
    this.extraStyles.defaultRuleString =
      " { fill: " + this.extraStyles.normalColor
      + "; rx: " + this.extraStyles.normalSize
      + "; ry: " + this.extraStyles.normalSize
      + "; }";

    // Begin drawing
    this.poincbox.unsuspendUpdate();
  };

  PoincareClickerController.prototype.handleTouch = function(b, lb, ptExists, thePoint) {
    if (!ptExists || !(thePoint.visProp.visible)) {
      if (insideSep(this.energy, b, lb)) {
        // console.log("("+b+","+lb+")=>");
        // console.log(initConds(this.energy, b, lb));
        var newPoincPoints = reapPoincarePoints(this.energy, this.npt, b, lb, this.deltaT, this.maxSteps);
        // console.log(newPoincPoints);

        var groupId = this.groupCounter;
        this.groupCounter++;
        var groupCSSClass = "pointGroup"+groupId;

        this.styleSheet.insertRule("."+groupCSSClass+ this.extraStyles.defaultRuleString);
        // Start with default color
        this.groupColMap[groupId] = 0;

        var overHandler = makePointOverHandler(this, groupId, groupCSSClass);
        var outHandler = makePointOutHandler(this, groupId, groupCSSClass);

        var newPointGroup = new Array(newPoincPoints.length);

        this.poincbox.suspendUpdate();
        for (var i = 0; i < newPoincPoints.length; i++) {
          newPointGroup[i] = this.poincbox.create('point', newPoincPoints[i],
                                                  this.basePointStyle);
          // newPointGroup[i].hasPoint = function(){return false; };
          // newPointGroup[i].rendNode.addEventListener('mouseenter', overHandler);
          // newPointGroup[i].rendNode.addEventListener('touchstart', overHandler);
          // newPointGroup[i].rendNode.addEventListener('mouseleave', outHandler);
          // newPointGroup[i].rendNode.addEventListener('touchend', outHandler);
          newPointGroup[i].on('over', overHandler);
          newPointGroup[i].on('out', outHandler);
          newPointGroup[i].groupId = groupId;
          newPointGroup[i].rendNode.classList.add(groupCSSClass);
        };
        this.poincbox.unsuspendUpdate();

        this.pointGroupList.push(newPointGroup);
      }
    } else {
      // point exists and is visible --- rotate color
      this.rotateColor(thePoint.groupId, thePoint.rendNode.classList[0]);

    };

    // QUESTION: Should we clear the 'Redo' stack?

    this.updateButtonAbility();

  };

  PoincareClickerController.prototype.rotateColor = function(groupId, groupCSSClass) {

    this.groupColMap[groupId] = (this.groupColMap[groupId] + 1) % this.colorPal.length;

    var selectorText = "." + groupCSSClass;
    var rules = this.styleSheet.cssRules;

    for(var i=0; i<rules.length; i++) {
      // Find the correct rule in the stylesheet
      if(rules[i].selectorText == selectorText) {
        rules[i].style['fill'] = this.colorPal[this.groupColMap[groupId]];
      };
    };

  };

  PoincareClickerController.prototype.handleKey = function(event) {
    if (event.repeat) // Ignore repeated keys
      return;
    switch (event.key.toLowerCase()) {
    case "z":
      this.popLastOrbit();
      break;
    case "x":
      this.restoreOrbit();
      break;
    case "a":
    case "m":
      this.morePointsFromLast();
      break;
    case "c":
      this.clearPoints();
      break;
    case "0":
      this.initialZoom();
      break;
    case "?":
      keyHelp();
      break;
    };
  };

  function keyHelp() {
    var intro = introJs();

    intro.setOptions({
      showStepNumbers: false,
      // showButtons: false,
      showBullets: false,
      overlayOpacity: 0.0,
      steps: [
        {
          element: "#poincbox",
          intro: "The following keyboard shortcuts are available:"+
            "<table><thead><tr><th>Key</th><th>Function</th></tr></thead><tbody>"+
            "<tr><td>?</td><td>Show this help</td></tr>"+
            "<tr><td>z</td><td>Undo add last orbit</td></tr>"+
            "<tr><td>x</td><td>Redo add last orbit</td></tr>"+
            "<tr><td>a or m</td><td>Extend last orbit</td></tr>"+
            "<tr><td>c</td><td>Clear board</td></tr>"+
            "<tr><td>0</td><td>Zoom 100%</td></tr>"+
            "</tbody></table>"
            ,
        }
      ]
    });

    intro.start();
  };

  PoincareClickerController.prototype.clearPoints = function() {

    if (this.poincbox.objectsList.length < this.clearPointThresh) {
      this._clearPointsSanely();
    } else { // JSXGraph is very slow at clearing a huge number of points
      this._clearPointsViaFreeBoard();
    };

    this.pointGroupList = new Array();
    this.undonePointGroupList = new Array();

    this.updateButtonAbility();
  };

  PoincareClickerController.prototype._clearPointsSanely = function() {

    this.poincbox.suspendUpdate();

    // Have to remove in reverse order to hack around JXG's
    // otherwise O(n^2) behaviour
    for (var i = this.poincbox.objectsList.length - 1;
         i >=0; i--) {
      // All the points on the board that we've added have the
      // property groupId
      if ((this.poincbox.objectsList[i].type == JXG.OBJECT_TYPE_POINT) &&
          this.poincbox.objectsList[i].hasOwnProperty('groupId'))
        this.poincbox.removeObject(this.poincbox.objectsList[i].id);
    };

    this.poincbox.unsuspendUpdate();
  };

  PoincareClickerController.prototype._clearPointsViaFreeBoard = function() {
    var bbox = this.poincbox.getBoundingBox();

    var poincboxName = this.poincbox.container;
    JXG.JSXGraph.freeBoard(this.poincbox);
    this.setupPoinc(poincboxName);
    this.poincbox.setBoundingBox(bbox);
  };

  PoincareClickerController.prototype.morePointsFromLast = function() {
    if (this.pointGroupList.length > 0) {
      var lastGroup = this.pointGroupList[this.pointGroupList.length - 1];
      if (lastGroup.length > 0) {
        var lastPoint = lastGroup[lastGroup.length - 1];
        var b = lastPoint.coords.usrCoords[1],
            lb = lastPoint.coords.usrCoords[2];

        var groupId = lastPoint.groupId;
        var groupCSSClass = "pointGroup"+groupId;
        var overHandler = makePointOverHandler(this, groupCSSClass);
        var outHandler = makePointOutHandler(this, groupCSSClass);

        var newPoincPoints = reapPoincarePoints(this.energy, this.npt, b, lb, this.deltaT, this.maxSteps);

        this.poincbox.suspendUpdate();
        for (var i = 0; i < newPoincPoints.length; i++) {
          var p = this.poincbox.create('point', newPoincPoints[i],
                                       this.basePointStyle);
          // p.hasPoint = function(){return false; };
          // p.rendNode.addEventListener('mouseenter', overHandler);
          // p.rendNode.addEventListener('touchstart', overHandler);
          // p.rendNode.addEventListener('mouseleave', outHandler);
          // p.rendNode.addEventListener('touchend', outHandler);
          p.on('over', overHandler);
          p.on('out', outHandler);
          p.groupId = groupId;
          p.rendNode.classList.add(groupCSSClass);
          // QUESTION: Is it better to resize the array instead of
          // constantly pushing?
          lastGroup.push(p);
        };
        this.poincbox.unsuspendUpdate();

      }
    }
  };

  PoincareClickerController.prototype.popLastOrbit = function() {
    if (this.pointGroupList.length == 0) return;

    var lastGroup = this.pointGroupList.pop();

    this.poincbox.suspendUpdate();
    for (var i = 0; i < lastGroup.length; i++) {
      lastGroup[i].hideElement();
    };

    this.undonePointGroupList.push(lastGroup);

    this.poincbox.unsuspendUpdate();

    this.updateButtonAbility();

  };

  PoincareClickerController.prototype.restoreOrbit = function() {
    if (this.undonePointGroupList.length == 0) return;

    var savedPoints = this.undonePointGroupList.pop();

    this.poincbox.suspendUpdate();
    for (var i = 0; i < savedPoints.length; i++) {
      savedPoints[i].showElement();
    };

    this.pointGroupList.push(savedPoints);

    this.poincbox.unsuspendUpdate();

    this.updateButtonAbility();
  };

  PoincareClickerController.prototype.bbox = function() {
    return [ -this.padFactor*bMax(this.energy), this.padFactor*lbMax(this.energy),
             this.padFactor*bMax(this.energy), -this.padFactor*lbMax(this.energy) ];
  };

  PoincareClickerController.prototype.initialZoom = function() {
    this.poincbox.suspendUpdate();
    this.poincbox.setBoundingBox(this.bbox(), false);
    this.poincbox.unsuspendUpdate();
    // This seems to be necessary to update what's shown on the axes
    this.poincbox.fullUpdate();

    this.isZoom100 = true;
    this.updateButtonAbility();
  }

  PoincareClickerController.prototype.updateButtonAbility = function() {
    this.undoButton.disabled = !(this.pointGroupList.length > 0);
    this.redoButton.disabled = !(this.undonePointGroupList.length > 0);
    this.zoomButton.disabled = this.isZoom100;
  };

  PoincareClickerController.prototype.setenergy = function(energy) {
    this.energy = energy;

    this.clearPoints();

    this.initialZoom();

    var oldnpt = this.npt;
    this.npt = this.nSlidePts;
    this.handleTouch(this.slideb, this.slidelb, false, null);
    this.npt = oldnpt;
  };

  PoincareClickerController.prototype.setnpt = function(npt) {
    this.npt = npt;
  };

  /* Add to global namespace */
  root['PoincareClickerController'] = PoincareClickerController;

  ////////////////////////////////////////////////////////////
  // For the torus demo

  // Constructor
  function TorusDemoController(threeboxName, ctrlsboxName) {

    if (!(this instanceof TorusDemoController)) {
      return new TorusDemoController(threeboxName, ctrlsboxName);
    }

    this.setupCtrls(ctrlsboxName);
    this.setupDemoGeom(threeboxName);
    this.setupEventListeners(threeboxName);

  };

  // Controller prototype
  TorusDemoController.prototype = {
    /* UI objects for the controls */
    ctrlsbox: {},
    freqslider: {},
    nptslider: {},

    /* threejs objects */
    three: {},
    trajObj: {},
    sectPointsObj: {},

    /* Variables */
    majorRad: 2.,
    minorRad: 0.6,
    trajBumpOutFac: 1.03,
    freqRatio: 1.618,
    deltaPhi: 0.04,
    nPhi: 48,
    nNewPhi: 8,
    nSectPoints: 14,

    tickTime: 35, // milliseconds
    timeoutID: {},
    runningState: false,

    /* Storage of trajectory points and section points */
    curPhi1: [0.],
    curPhi2: [0.],
    sectPoints: new Array(),

    /* Public member functions */
    setupCtrls: {},
    setupDemoGeom: {},
    setupEventListeners: {},
    advanceTraj: {},
    anglesTo3d: {},
    sectPointsToAdd: {},
    runTick: {},
    toggleRunning: {},

  };

  TorusDemoController.prototype.setupCtrls = function(ctrlsboxName) {
    this.ctrlsbox = JXG.JSXGraph.initBoard(ctrlsboxName,
                                 {boundingbox:[0.,1.,1.,0.],
                                  axis:false,
                                  pan: {enabled: false},
                                  showNavigation: false,
                                  showCopyright:  false});
    this.ctrlsbox.suspendUpdate();

    this.freqslider = this.ctrlsbox.create(
      'slider',
      [[0.05,.66],[0.63,.66],
       [0.5, 1.618, 2.]],
      {name: 'Frequency ratio', precision:3});

    this.nptslider = this.ctrlsbox.create(
      'slider',
      [[0.05,.33],[0.63,.33],
       [1,this.nSectPoints,100]],
      {name: '# pts on section', snapWidth:1, precision:0});

    this.freqslider.on('drag',
                       (function (o) {
                         return function() {
                           o.freqRatio = o.freqslider.Value();
                           o.sectPoints.splice(0, o.sectPoints.length-1);
                         };
                       })(this));
    this.nptslider.on('drag',
                      (function (o) {
                         return function() {
                           o.nSectPoints = o.nptslider.Value();
                         };
                       })(this));

    this.ctrlsbox.unsuspendUpdate();
  };

  TorusDemoController.prototype.setupDemoGeom = function(threeboxName) {
    // Bootstrap core + controls plugin into element
    var three = THREE.Bootstrap({
      element: document.querySelector('#' + threeboxName),
      plugins: ['core', 'controls', 'cursor'],
      controls: {
        klass: THREE.OrbitControls
      },
    });

    this.three = three;

    three.scene.add( new THREE.AmbientLight( 0xf0f0f0 ) );

    var directionalLight = new THREE.DirectionalLight( 0xf0f0f0, 0.5 );
    three.scene.add( directionalLight );

    var helper = new THREE.GridHelper( 60, 5 );
    helper.position.y = - 6;
    helper.material.opacity = 0.25;
    helper.material.transparent = true;
    three.scene.add( helper );

    three.renderer.setClearColor( 0xf0f0f0 );

    // Torus
    var geometry = new THREE.TorusGeometry( this.majorRad, this.minorRad, 16, 50 );
    geometry.rotateX(Math.PI * 0.5);
    var material = new THREE.MeshStandardMaterial( {
      transparent: true, opacity: 0.35,
      emissive: 0x0, roughness: 0.2, metalness: 0.2,
      side: THREE.DoubleSide,
      wireframe: false,
      depthTest: false, depthWrite: false,
      color: 0x21ce70,  } );
    var torus = new THREE.Mesh( geometry, material );
    three.scene.add( torus );

    // Cutting section
    var sectSize = 2. * this.majorRad;

    var plane = new THREE.PlaneBufferGeometry( sectSize, sectSize );
    var planeMat = new THREE.MeshStandardMaterial( {
      transparent: true, opacity: 0.25,
      emissive: 0x0, roughness: 0.2, metalness: 0.2,
      wireframe: false,
      depthTest: false, depthWrite: false,
      side: THREE.DoubleSide,
      color: 0x0000ff,  } );

    var sect = new THREE.Mesh( plane, planeMat );

    sect.position.x = this.majorRad;
    three.scene.add( sect );

    // Alter controls
    three.controls.rotateSpeed = 0.5;

    // Place camera
    three.camera.position.set(3, 2.25, 5.);

    // For the trajectory object
    this.trajObj = new THREE.Line( new THREE.Geometry(),
                                   new THREE.LineBasicMaterial( { color : 0x000000,
                                                                  linewidth : 2 } ));

    three.scene.add( this.trajObj );

    // For the poinc sect points
    var sectPointGeom = new THREE.Geometry();
    var minor_r = this.minorRad * this.trajBumpOutFac;
    sectPointGeom.vertices.push( new THREE.Vector3( this.majorRad + minor_r, 0., 0.) );

    var pointMat = new THREE.PointsMaterial( { size: 5,
                                               color: 0x880088,
                                               sizeAttenuation: false,
                                               alphaTest: 0.5,
                                               transparent: true } );

    this.sectPointsObj = new THREE.Points( sectPointGeom, pointMat );
    three.scene.add( this.sectPointsObj );

    this.advanceTraj();

  };

  TorusDemoController.prototype.setupEventListeners = function(threeboxName) {

    var threeBoxEl = document.getElementById(threeboxName);

    window.addEventListener('scroll',
                            (function(o, el) {
                              return function(e) {

                                var isVis = isPartiallyScrolledIntoView(el);

                                if (o.runningState) {
                                  if (!isVis) {
                                    o.toggleRunning();
                                  }
                                } else {
                                  // not yet running
                                  if (isVis) {
                                    o.toggleRunning();
                                  }
                                }
                              }
                            })(this, threeBoxEl));
  };

  TorusDemoController.prototype.anglesTo3d = function(phi1, phi2) {
    var n_phi = phi1.length;
    var threeVectors = new Array(n_phi);

    var minor_r = this.minorRad * this.trajBumpOutFac;

    for (var i=0; i<n_phi; i++) {
      var rho = this.majorRad + minor_r * Math.cos(phi2[i]);

      // graphics people... argh
      threeVectors[i] = new THREE.Vector3(
        rho * Math.cos(phi1[i]),
        minor_r * Math.sin(phi2[i]),
        rho * Math.sin(phi1[i]));
    };

    return threeVectors;
  };

  // Because mods should always be in the range 0 <= x < m
  function fmod(x, m)
  {
    return ((x<0)? Math.abs(m) : 0) + (x % m);
  };

  TorusDemoController.prototype.sectPointsToAdd = function(vects) {
    var n_vect = vects.length;

    var sectPoints = new Array();

    for (var i=1; i<n_vect; i++) {
      var v0 = vects[i-1];
      var v1 = vects[i];
      if ((v1.x > 0) && (v0.z < 0) && (v1.z >= 0)) {
        // linearly interpolate between v0 and v1, using the z coord as the indep var
        var newPoint = new THREE.Vector3(
          linearInterp(0., v0.z, v0.x, v1.z, v1.x),
          linearInterp(0., v0.z, v0.y, v1.z, v1.y),
          0.
        );
        sectPoints.push(newPoint);
      };
    };

    return sectPoints;

  };

  TorusDemoController.prototype.advanceTraj = function() {
    var phi1 = new Array(this.nNewPhi);
    var phi2 = new Array(this.nNewPhi);

    var phi10 = fmod(this.curPhi1[ this.curPhi1.length - 1 ],
                     2. * Math.PI);
    var phi20 = fmod(this.curPhi2[ this.curPhi2.length - 1 ],
                     2. * Math.PI);

    for (var i=0; i<this.nNewPhi; i++) {
      phi1[i] = phi10 + i * this.deltaPhi;
      phi2[i] = phi20 + i * this.deltaPhi * this.freqRatio;
    };

    var newPhi1 = this.curPhi1.concat(phi1);
    var newPhi2 = this.curPhi2.concat(phi2);

    // Drop enough to keep the desired length
    newPhi1.splice(0, newPhi1.length - this.nPhi);
    newPhi2.splice(0, newPhi2.length - this.nPhi);

    this.curPhi1 = newPhi1;
    this.curPhi2 = newPhi2;

    var threeVectors = this.anglesTo3d(newPhi1, newPhi2);


    var newGeom = new THREE.Geometry();
    newGeom.vertices = threeVectors;

    this.trajObj.geometry.dispose();
    this.trajObj.geometry = newGeom;

    // Find possible points to add to section
    // They can only come from the last nNewPhi points, which are the new ones
    var lastVects = threeVectors.slice(this.nPhi - this.nNewPhi);
    var newSectPoints = this.sectPointsToAdd(lastVects);

    newSectPoints = this.sectPoints.concat(newSectPoints);
    newSectPoints.splice(0, newSectPoints.length - this.nSectPoints);

    this.sectPoints = newSectPoints;

    var newSectGeom = new THREE.Geometry();
    newSectGeom.vertices = newSectPoints;

    if (newSectPoints.length > 0) {
      this.sectPointsObj.material.visible = true;
      this.sectPointsObj.geometry.dispose();
      this.sectPointsObj.geometry = newSectGeom;
    } else {
      this.sectPointsObj.material.visible = false;
    };

  };

  TorusDemoController.prototype.toggleRunning = function() {
    if (this.runningState) {
      window.clearTimeout(this.timeoutID);
      this.runningState = false;
    } else {
      this.runningState = true;
      // Will set its own timeout
      this.runTick();
    };
  };

  TorusDemoController.prototype.runTick = function() {
    this.advanceTraj();

    if (this.runningState) {
      // OMGauss I hate javascript...
      this.timeoutID = setTimeout((function(o) {
                                     return function() { o.runTick() };
                                  })(this),
                                  this.tickTime);
    } else {
      // nothing to do
    };
  };

  /* Add to global namespace */
  root['TorusDemoController'] = TorusDemoController;

  ////////////////////////////////////////////////////////////

  /* For tour using intro.js */
  root['startIntro'] = function(){
    var intro = introJs();

    intro.setOptions({
      disableInteraction: false,
      overlayOpacity: 0.0,
      steps: [
        {
          element: '#poincbox',
          intro: "This is where the Poincaré section lives.  Clicking inside the box will add points from a phase space orbit."
        },
        {
          element: '#ctrlsbox',
          intro: "Here is where you set the energy of the system, and how many points are added on each click.  Changing the energy will clear all points.",
        },
        {
          element: '#buttonbox',
          intro: "Use 'Extend last orbit' to add n more points to the most recent orbit.",
        },
        {
          element: '#poincbox',
          intro: "When you hover over a point, all points from that orbit will be highlighted.  Clicking on an orbit will assign a new color, so you can keep track of different orbits.  You can Shift-drag to pan, and Ctrl-drag a region to zoom in for more detail.",
        },
      ]
    });

    intro.start();
  };


})(this);
