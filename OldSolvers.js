//3456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_
// (JT: why the numbers? counts columns, helps me keep 80-char-wide listings)

/*=====================
  VBOboxes.js library: 
  ===================== 
One 'VBObox' object contains all we need for WebGL/OpenGL to render on-screen a 
		the shapes made from vertices stored in one Vertex Buffer Object (VBO), 
		as drawn by one 'shader program' that runs on your computer's Graphical  
		Processing Unit(GPU).  
The 'shader program' consists of a Vertex Shader and a Fragment Shader written 
		in GLSL, compiled and linked and ready to execute as a Single-Instruction, 
		Multiple-Data (SIMD) parallel program executed simultaneously by multiple 
		'shader units' on the GPU.  The GPU runs one 'instance' of the Vertex 
		Shader for each vertex in every shape, and one 'instance' of the Fragment 
		Shader for every on-screen pixel covered by any part of any drawing 
		primitive defined by those vertices.
The 'VBO' consists of a 'buffer object' (a memory block reserved in the GPU),
		accessed by the shader program through its 'attribute' and 'uniform' 
		variables.  Each VBObox object stores its own 'uniform' values as vars in
		JavaScript; its 'adjust()'	function computes newly-updated values for these
		vars and then transfers them to the GPU for use.
EVENTUALLY you should replace 'cuon-matrix' with the free, open-source
   'glmatrix.js' library for vectors, matrices & quaternions: Google it!
		This vector/matrix library is more complete, more widely-used, and runs
		faster than our textbook's 'cuon-matrix' library.  
		--------------------------------------------------------------
		I recommend you use glMatrix.js instead of cuon-matrix-quat03.js
		--------------------------------------------------------------
		for all future WebGL programs. You can CONVERT existing cuon-matrix-based
		programs to glmatrix.js in a very gradual, sensible, testable way:
		--add the glmatrix.js library to an existing cuon-matrix-based program;
			(but don't call any of its functions yet)
		--comment out the glmatrix.js parts (if any) that cause conflicts or in	
			any way disrupt the operation of your program.
		--make just one small local change in your program; find a small, simple,
			easy-to-test portion of your program where you can replace a 
			cuon-matrix object or function call with a glmatrix function call.
			Test; make sure it works.
		--Save a copy of this new program as your latest numbered version. Repeat
			the previous step: go on to the next small local change in your program
			and make another replacement of cuon-matrix use with glmatrix use. 
			Test it; make sure it works; save this as your next numbered version.
		--Continue this process until your program no longer uses any cuon-matrix
			library features at all, and no part of glmatrix is commented out.
			Remove cuon-matrix from your library, and now use only glmatrix.

	-------------------------------------------------------
	A MESSY SET OF CUSTOMIZED OBJECTS--NOT REALLY A 'CLASS'
	-------------------------------------------------------
As each 'VBObox' object will contain DIFFERENT GLSL shader programs, DIFFERENT 
		attributes for each vertex, DIFFERENT numbers of vertices in VBOs, and 
		DIFFERENT uniforms, I don't see any easy way to use the exact same object 
		constructors and prototypes for all VBObox objects.  Every additional VBObox 
		objects may vary substantially, so I recommend that you copy and re-name an 
		existing VBObox prototype object, and modify as needed, as shown here. 
		(e.g. to make the VBObox2 object, copy the VBObox1 constructor and 
		all its prototype functions, then modify their contents for VBObox2 
		activities.)
Note that you don't really need multiple 'VBObox' objects for any simple, 
		beginner-level WebGL/OpenGL programs: if all vertices contain exactly 
		the same attributes (e.g. position, color, surface normal), and use 
		the same shader program (e.g. same Vertex Shader and Fragment Shader), 
		then our textbook's simple 'example code' will suffice.  
		But that's rare -- most genuinely useful WebGL/OpenGL programs need 
		different sets of vertices with  different sets of attributes rendered 
		by different shader programs.  
		*** A customized VBObox object for each VBO/shader pair will help you
		remember and correctly implement ALL the WebGL/GLSL steps required for 
		a working multi-shader, multi-VBO program.
*/
// Written for EECS 351-2,	Intermediate Computer Graphics,
//							Northwestern Univ. EECS Dept., Jack Tumblin
// 2016.05.26 J. Tumblin-- Created; tested on 'TwoVBOs.html' starter code.
// 2017.02.20 J. Tumblin-- updated for EECS 351 use for Project C.
// 2018.04.09 J. Tumblin-- minor corrections/renaming for particle systems.
//=============================================================================
// Tabs set to 2

//=============================================================================
//=============================================================================
function VBObox0() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox0' object  that holds all data and 
// fcns needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate set of shaders.

    this.VERT_SRC =	//--------------------- VERTEX SHADER source code
        'precision highp float;\n' +				// req'd in OpenGL ES if we use 'float'
        //
        'uniform mat4 u_MvpMatrix;\n' +
        'attribute vec4 a_Position;\n' +
        'attribute vec3 a_Color;\n' +
        //'attribute float a_PtSize; \n' +
        'varying vec3 v_Colr2;\n' +
        //
        'void main() {\n' +
        '  gl_PointSize = 20.0;\n' +
        '  gl_Position = u_MvpMatrix *  a_Position;\n' +
        '	 v_Colr2 = a_Color;\n' +
        ' }\n';

    this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code
        'precision mediump float;\n' +
        'varying vec3 v_Colr2;\n' +
        'void main() {\n' +
        // '  float dist = distance(gl_PointCoord, vec2(0.5, 0.5)); \n' +
        // '  if(dist < 0.5) {\n' +
        '  	gl_FragColor = vec4(v_Colr2.rgb, 1.0);\n' +
        //	'  	gl_FragColor = vec4((1.0-2.0*dist)*v_Colr2.rgb, 1.0);\n' +
        //'    } else {discard;};' +
        '}\n';


    this.vboVerts = gndVerts.length / floatsPerVertex;							// # of vertices held in 'vboContents' array;
    this.vboLoc;										// Vertex Buffer Object location# on the GPU
    this.FSIZE = gndVerts.BYTES_PER_ELEMENT;
    // bytes req'd for 1 array element;
    // (why? used to compute stride and offset
    // in bytes for vertexAttribPointer() calls)
    this.shaderLoc;									// Shader-program location # on the GPU, made
    // by compile/link of VERT_SRC and FRAG_SRC.
    //-------------------- Attribute locations in our shaders
    this.a_PositionLoc;							// GPU location: shader 'a_Position' attribute
    this.a_ColorLoc;								// GPU location: shader 'a_Color' attribute
    this.a_PtSizeLoc;								// GPU location: shader 'a_PtSize' attribute
    //-------------------- Uniform locations &values in our shaders
    this.MvpMatrix = new Matrix4();	// Transforms CVV axes to model axes.
    this.u_MvpMatrixLoc;						// GPU location for u_ModelMat uniform

}

VBObox0.prototype.init = function (myGL) {//=============================================================================
// Create, compile, link this VBObox object's shaders to an executable 'program'
// ready for use in the GPU.  Create and fill a Float32Array that holds all VBO 
// vertices' values; create a new VBO on the GPU and fill it with those values. 
// Find the GPU location of	all our shaders' attribute- and uniform-variables; 
// assign the correct portions of VBO contents as the data source for each 
// attribute, and transfer current values to the GPU for each uniform variable.
// (usually called only once, within main()) 

    this.shaderLoc = createProgram(myGL, this.VERT_SRC, this.FRAG_SRC);
    if (!this.shaderLoc) {
        console.log(this.constructor.name +
            '.init() failed to create executable Shaders on the GPU. Bye!');
        return;
    }
    this.vboLoc = myGL.createBuffer();
    if (!this.vboLoc) {
        console.log(this.constructor.name +
            '.init() failed to create VBO in GPU. Bye!');
        return;
    }
    // Specify the purpose of our newly-created VBO.  Your choices are:
    //	== "gl.ARRAY_BUFFER" : the VBO holds vertices, each made of attributes
    // (positions, colors, normals, etc), or
    //	== "gl.ELEMENT_ARRAY_BUFFER" : the VBO holds indices only; integer values
    // that each select one vertex from a vertex array stored in another VBO.
    myGL.bindBuffer(myGL.ARRAY_BUFFER,	// GLenum 'target' for this GPU buffer
        this.vboLoc);				// the ID# the GPU uses for this buffer.

    // Transfer data from our JavaScript Float32Array object to the just-bound VBO.
    //  (Recall gl.bufferData() changes GPU's memory allocation: use
    //		gl.bufferSubData() to modify buffer contents without changing its size)
    //	The 'hint' helps GPU allocate its shared memory for best speed & efficiency
    //	(see OpenGL ES specification for more info).  Your choices are:
    //		--STATIC_DRAW is for vertex buffers rendered many times, but whose
    //				contents rarely or never change.
    //		--DYNAMIC_DRAW is for vertex buffers rendered many times, but whose
    //				contents may change often as our program runs.
    //		--STREAM_DRAW is for vertex buffers that are rendered a small number of
    // 			times and then discarded; for rapidly supplied & consumed VBOs.
    myGL.bufferData(gl.ARRAY_BUFFER, 			// GLenum target(same as 'bindBuffer()')
        gndVerts, 		// JavaScript Float32Array
        gl.STATIC_DRAW);			// Usage hint.

// Find & Set All Attributes:------------------------------
    // a) Get the GPU location for each attribute var used in our shaders:
    this.a_PositionLoc = gl.getAttribLocation(this.shaderLoc, 'a_Position');
    if (this.a_PositionLoc < 0) {
        console.log(this.constructor.name +
            '.init() Failed to get GPU location of attribute a_Position');
        return -1;	// error exit.
    }
    this.a_ColorLoc = myGL.getAttribLocation(this.shaderLoc, 'a_Color');
    if (this.a_ColorLoc < 0) {
        console.log(this.constructor.name +
            '.init() failed to get the GPU location of attribute a_Color');
        return -1;	// error exit.
    }
    // NEW! a_PtSize' attribute values are stored only in VBO2, not VBO1:
    // this.a_PtSizeLoc = gl.getAttribLocation(this.shaderLoc, 'a_PtSize');
    // if(this.a_PtSizeLoc < 0) {
    //   console.log(this.constructor.name +
//	    					'.init() failed to get the GPU location of attribute a_PtSize');
//	  return -1;	// error exit.
    // }
    // b) Next, set up GPU to fill these attribute vars in our shader with
    // values pulled from the currently-bound VBO (see 'gl.bindBuffer()).
    // 	Here's how to use the almost-identical OpenGL version of this function:
    //		http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribPointer.xml )
    myGL.vertexAttribPointer(
        this.a_PositionLoc,//index == ID# for the attribute var in GLSL shader pgm;
        4,						// size == how many dimensions for this attribute: 1,2,3 or 4?
        myGL.FLOAT,		// type == what data type did we use for those numbers?
        false,				// isNormalized == are these fixed-point values that we need
        //									normalize before use? true or false
        7 * this.FSIZE,	// Stride == #bytes we must skip in the VBO to move from one
        // of our stored attributes to the next.  This is usually the
        // number of bytes used to store one complete vertex.  If set
        // to zero, the GPU gets attribute values sequentially from
        // VBO, starting at 'Offset'.
        // (Vertex size in bytes: 4 floats for pos; 3 color; 1 PtSize)
        0);						// Offset == how many bytes from START of buffer to the first
    // value we will actually use?  (We start with position).
    myGL.vertexAttribPointer(this.a_ColorLoc, 3, myGL.FLOAT, false,
        7 * this.FSIZE, 			// stride for VBO2 (different from VBO1!)
        4 * this.FSIZE);			// offset: skip the 1st 4 floats.
//  myGL.vertexAttribPointer(this.a_PtSizeLoc, 1, myGL.FLOAT, false, 
//							8*this.FSIZE,		// stride for VBO2 (different from VBO1!) 
//							7*this.FSIZE);		// offset: skip the 1st 7 floats.

    // c) Enable this assignment of the attribute to its' VBO source:
    myGL.enableVertexAttribArray(this.a_PositionLoc);
    myGL.enableVertexAttribArray(this.a_ColorLoc);

// Find All Uniforms:--------------------------------
//Get GPU storage location for each uniform var used in our shader programs: 
    this.u_MvpMatrixLoc = gl.getUniformLocation(this.shaderLoc, 'u_MvpMatrix');
    if (!this.u_MvpMatrixLoc) {
        console.log(this.constructor.name +
            '.init() failed to get GPU location for u_MvpMatrix uniform');
        return;
    }

}

VBObox0.prototype.adjust = function (myGL) {
//=============================================================================
// Update the GPU to newer, current values we now store for 'uniform' vars on 
// the GPU; and (if needed) update each attribute's stride and offset in VBO.
    myGL.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    /*																	// -and-linked executable shader program.
      // Adjust values for our uniforms: -----------------------
    this.ModelMatrix.setRotate(-g_currentAngle, 0, 0, 1);	// -spin drawing axes,
    this.ModelMatrix.translate(0.35, -0.15, 0);						// then translate them.
    //  Transfer new uniforms' values to the GPU:-------------*/
    // Send  new 'ModelMat' values to the GPU's 'u_ModelMat1' uniform:
    myGL.viewport(0,											 				// Viewport lower-left corner
        0, 			// location(in pixels)
        myGL.drawingBufferWidth, 					// viewport width,
        myGL.drawingBufferHeight);			// viewport height in pixels.


    this.vpAspect = ((myGL.drawingBufferWidth) /			// On-screen aspect ratio for
        myGL.drawingBufferHeight);		// this camera: width/height.

    this.MvpMatrix.setIdentity();	// rotate drawing axes,
    this.MvpMatrix.perspective(45, this.vpAspect, 1, +200);							// then translate them.
    this.MvpMatrix.lookAt(xCamPos, yCamPos, zCamPos, xCamPos + Math.cos(Math.PI * (currentAngle / 180)),
        yCamPos + Math.sin(Math.PI * (currentAngle / 180)), zCamPos + zlook, 0.0, 0, 1);   // UP vector.

    myGL.uniformMatrix4fv(this.u_MvpMatrixLoc,	// GPU location of the uniform
        false, 										// use matrix transpose instead?
        this.MvpMatrix.elements);	// send data from Javascript.
}
/*
VBObox0.prototype.PartSys_render = function(myGL) {

//=============================================================================
// Send commands to GPU to select and render current VBObox contents.  	
 if(isClear == 1) gl.clear(gl.COLOR_BUFFER_BIT);
  myGL.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
  																	// -and-linked executable shader program.
//------CAREFUL! RE-BIND YOUR VBO AND RE-ASSIGN SHADER ATTRIBUTES!-------------
//		Each call to useProgram() reconfigures the GPU's processors & data paths 
// for efficient SIMD execution of the newly-selected shader program. While the 
// 'old' shader program's attributes and uniforms remain at their same memory 
// locations, starting the new shader program invalidates the old data paths 
// that connected these attributes to the VBOs in memory that supplied their 
// values. When we call useProgram() to return to our 'old' shader program, we 
// must re-establish those data-paths between shader attributes and VBOs, even 
// if those attributes, VBOs, and locations have not changed!
//		Thus after each useProgram() call, we must:
// a)--call bindBuffer() again to re-bind each VBO that our shader will use, &
// b)--call vertexAttribPointer() again for each attribute in our new shader
//		program, to re-connect the data-path(s) from bound VBO(s) to attribute(s):
// c)--call enableVertexAttribArray() to enable use of those data paths.
//----------------------------------------------------
	// a) Re-set the GPU's currently 'bound' vbo buffer;
	myGL.bindBuffer(myGL.ARRAY_BUFFER,	// GLenum 'target' for this GPU buffer 
										this.vboLoc);			// the ID# the GPU uses for this buffer.
	// (Here's how to use the almost-identical OpenGL version of this function:
	//		http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribPointer.xml )
  // b) Re-connect data paths from VBO to each shader attribute:
  myGL.vertexAttribPointer(this.a_PosLoc, 3, myGL.FLOAT, false, 
  													PART_MAXVAR*this.FSIZE, 0);		// stride, offset
 myGL.vertexAttribPointer(this.a_ColrLoc, 3, myGL.FLOAT, false, 
  												PART_MAXVAR*this.FSIZE, PART_R*this.FSIZE); // stride, offset
 myGL.vertexAttribPointer(this.a_diamID, 1, myGL.FLOAT, false, 
  												PART_MAXVAR*this.FSIZE, PART_DIAM*this.FSIZE); // stride, offset
  // c) enable the newly-re-assigned attributes:
  myGL.enableVertexAttribArray(this.a_PosLoc);
	myGL.enableVertexAttribArray(this.a_ColrLoc);
	myGL.enableVertexAttribArray(this.a_diamID);

  // ----------------------------Draw the contents of the currently-bound VBO:

	
	myGL.bufferSubData(myGL.ARRAY_BUFFER, 0, s0);
  myGL.drawArrays(myGL.POINTS, 0, partCount);
}

*/

VBObox0.prototype.draw = function (myGL) {
//=============================================================================
// Send commands to GPU to select and render current VBObox contents.
    if (isClear == 1) gl.clear(gl.COLOR_BUFFER_BIT);

    myGL.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    // -and-linked executable shader program.
//------CAREFUL! RE-BIND YOUR VBO AND RE-ASSIGN SHADER ATTRIBUTES!-------------
//		Each call to useProgram() reconfigures the GPU's processors & data paths 
// for efficient SIMD execution of the newly-selected shader program. While the 
// 'old' shader program's attributes and uniforms remain at their same memory 
// locations, starting the new shader program invalidates the old data paths 
// that connected these attributes to the VBOs in memory that supplied their 
// values. When we call useProgram() to return to our 'old' shader program, we 
// must re-establish those data-paths between shader attributes and VBOs, even 
// if those attributes, VBOs, and locations have not changed!
//		Thus after each useProgram() call, we must:
// a)--call bindBuffer() again to re-bind each VBO that our shader will use, &
// b)--call vertexAttribPointer() again for each attribute in our new shader
//		program, to re-connect the data-path(s) from bound VBO(s) to attribute(s):
// c)--call enableVertexAttribArray() to enable use of those data paths.
//----------------------------------------------------
    // a) Re-set the GPU's currently 'bound' vbo buffer;
    myGL.bindBuffer(myGL.ARRAY_BUFFER,	// GLenum 'target' for this GPU buffer
        this.vboLoc);			// the ID# the GPU uses for this buffer.
    // (Here's how to use the almost-identical OpenGL version of this function:
    //		http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribPointer.xml )
    //b) Re-connect data paths from VBO to each shader attribute:
    myGL.vertexAttribPointer(this.a_PositionLoc, 4, myGL.FLOAT, false,
        7 * this.FSIZE, 0);							// stride, offset
    myGL.vertexAttribPointer(this.a_ColorLoc, 3, myGL.FLOAT, false,
        7 * this.FSIZE, 4 * this.FSIZE);	// stride, offset
//	myGL.vertexAttribPointer(this.a_PtSizeLoc,   1, myGL.FLOAT, false,
//														8*this.FSIZE, 7*this.FSIZE);	// stride, offset
    // c) Re-Enable use of the data path for each attribute:
    myGL.enableVertexAttribArray(this.a_PositionLoc);
    myGL.enableVertexAttribArray(this.a_ColorLoc);
    // ----------------------------Draw the contents of the currently-bound VBO:
    myGL.drawArrays(myGL.LINES, 		// select the drawing primitive to draw,
        0, 								// location of 1st vertex to draw;
        this.vboVerts);		// number of vertices to draw on-screen.
}
//  GROUND GRID


VBObox1.prototype.reset = function () {
    /*for (i = 0; i < partCount; i++){
                if(s0[i* PART_MAXVAR + PART_XVEL] > 0.0) s0[i* PART_MAXVAR + PART_XVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_XVEL] -= INIT_VEL;
                if(s0[i* PART_MAXVAR + PART_YVEL] > 0.0) s0[i* PART_MAXVAR + PART_YVEL] += INIT_VEL*.8; else s0[i* PART_MAXVAR + PART_YVEL] -= INIT_VEL*.8;
                if(s0[i* PART_MAXVAR + PART_ZVEL] > 0.0) s0[i* PART_MAXVAR + PART_ZVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_ZVEL] -= INIT_VEL;
                }*/
    for (var i = 0; i < partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (this.s0[pOff + PART_XVEL] > 0) {
            this.s0[pOff + PART_XVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else this.s0[pOff + PART_XVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (this.s0[pOff + PART_YVEL] > 0) {
            this.s0[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else this.s0[pOff + PART_YVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (this.s0[pOff + PART_ZVEL] > 0) {
            this.s0[pOff + PART_ZVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else this.s0[pOff + PART_ZVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
    }
}
//===================Mouse and Keyboard event-handling Callbacks================
//==============================================================================
VBObox0.prototype.myMouseDown = function (ev, gl, canvas) {
//==============================================================================
// Called when user PRESSES down any mouse button;
// 									(Which button?    console.log('ev.button='+ev.button);   )
// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)  

// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
    var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
    var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
    var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
//  console.log('myMouseDown(pixel coords): xp,yp=\t',xp,',\t',yp);

    // Convert to Canonical View Volume (CVV) coordinates too:
    var x = (xp - canvas.width / 2) / 		// move origin to center of canvas and
        (canvas.width / 2);			// normalize canvas to -1 <= x < +1,
    var y = (yp - canvas.height / 2) /		//										 -1 <= y < +1.
        (canvas.height / 2);
//	console.log('myMouseDown(CVV coords  ):  x, y=\t',x,',\t',y);

    isDrag = true;											// set our mouse-dragging flag
    xMclik = x;													// record where mouse-dragging began
    yMclik = y;
    if (x < 0) ANGLE3_STEP = 45;
    else ANGLE3_STEP = -45;
    document.getElementById('MouseResult1').innerHTML =
        'myMouseDown() at CVV coords x,y = ' + x + ', ' + y + '<br>';
};


VBObox0.prototype.myMouseMove = function (ev, gl, canvas) {
//==============================================================================
// Called when user MOVES the mouse with a button already pressed down.
// 									(Which button?   console.log('ev.button='+ev.button);    )
// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)  

    if (isDrag == false) return;				// IGNORE all mouse-moves except 'dragging'

    // Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
    var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
    var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
    var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
//  console.log('myMouseMove(pixel coords): xp,yp=\t',xp,',\t',yp);

    // Convert to Canonical View Volume (CVV) coordinates too:
    var x = (xp - canvas.width / 2) / 		// move origin to center of canvas and
        (canvas.width / 2);			// normalize canvas to -1 <= x < +1,
    var y = (yp - canvas.height / 2) /		//										 -1 <= y < +1.
        (canvas.height / 2);
//	console.log('myMouseMove(CVV coords  ):  x, y=\t',x,',\t',y);

    // find how far we dragged the mouse:
    xMdragTot += (x - xMclik);					// Accumulate change-in-mouse-position,&
    yMdragTot += (y - yMclik);
    xMclik = x;													// Make next drag-measurement from here.
    yMclik = y;
    //if((x-xMclik) > 0) ANGLE3_STEP = 60;
    //else ANGLE3_STEP = -60;
// (? why no 'document.getElementById() call here, as we did for myMouseDown()
// and myMouseUp()? Because the webpage doesn't get updated when we move the 
// mouse. Put the web-page updating command in the 'draw()' function instead)
};

VBObox0.prototype.myMouseUp = function (ev, gl, canvas) {
//==============================================================================
// Called when user RELEASES mouse button pressed previously.
// 									(Which button?   console.log('ev.button='+ev.button);    )
// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)  

// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
    ANGLE3_STEP = 0;
    var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
    var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
    var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
//  console.log('myMouseUp  (pixel coords): xp,yp=\t',xp,',\t',yp);

    // Convert to Canonical View Volume (CVV) coordinates too:
    var x = (xp - canvas.width / 2) / 		// move origin to center of canvas and
        (canvas.width / 2);			// normalize canvas to -1 <= x < +1,
    var y = (yp - canvas.height / 2) /		//										 -1 <= y < +1.
        (canvas.height / 2);
    console.log('myMouseUp  (CVV coords  ):  x, y=\t', x, ',\t', y);

    isDrag = false;											// CLEAR our mouse-dragging flag, and
    // accumulate any final bit of mouse-dragging we did:
    xMdragTot += (x - xMclik);
    yMdragTot += (y - yMclik);
    console.log('myMouseUp: xMdragTot,yMdragTot =', xMdragTot, ',\t', yMdragTot);
    // Put it on our webpage too...
    document.getElementById('MouseResult1').innerHTML =
        'myMouseUp(       ) at CVV coords x,y = ' + x + ', ' + y + '<br>';
};


VBObox0.prototype.myKeyDown = function (ev) {


    if (ev.keyCode == 68)  // D
        ANGLE3_STEP = -60;
    if (ev.keyCode == 65) // A
        ANGLE3_STEP = 60;
    if (ev.keyCode == 87) // W
        zStep = 10;
    if (ev.keyCode == 83)  // S
        zStep = -10;
    if (ev.keyCode == 73) //I
        zMoveStep = 1;
    if (ev.keyCode == 75)  //K
        zMoveStep = -1;

    if (ev.keyCode == 38) {  //up
        moveSpeed = .1;
        return;
    }
    if (ev.keyCode == 40) {  //Down
        moveSpeed = -.1;
        return;
    }
    if (ev.keyCode == 37) {  //Left
        moveSpeed2 = -.1;
        return;
    }
    if (ev.keyCode == 39) {  //Right
        moveSpeed2 = .1;
        return;
    }
    if (ev.keyCode == 89) { // Y
        springRun = 1;
    }
}

VBObox0.prototype.myKeyUp = function (ev) {

    if (ev.keyCode == 68)
        ANGLE3_STEP = 0;
    if (ev.keyCode == 65)
        ANGLE3_STEP = 0;
    if (ev.keyCode == 87)
        zStep = 0;
    if (ev.keyCode == 83)
        zStep = -0;
    if (ev.keyCode == 73)
        zMoveStep = 0;
    if (ev.keyCode == 75)
        zMoveStep = 0;
    if (ev.keyCode == 38) {
        moveSpeed = 0;
        return;
    }
    if (ev.keyCode == 40) {
        moveSpeed = 0;
        return;
    }
    if (ev.keyCode == 37) {
        moveSpeed2 = 0;
        return;
    }
    if (ev.keyCode == 39) {
        moveSpeed2 = 0;
        return;
    }
    if (ev.keyCode == 89) {
        springRun = 0;
    }
}

VBObox0.prototype.myKeyPress = function (ev, partSys) {

    myChar = String.fromCharCode(ev.keyCode);	//	convert code to character-string

    switch (myChar) {
        case '0':
            myRunMode = 0;			// RESET!
            break;
        case '1':
            myRunMode = 1;			// PAUSE!
            break;
        case '2':
            myRunMode = 2;			// STEP!
            break;
        case '3':							// RUN!
            myRunMode = 3;
            break;
        case 'b':							// Toggle floor-bounce constraint type:
        case 'B':
            if (g_bounce == 0) g_bounce = 1;
            else g_bounce = 0;
            break;
        case 'c':					// 'c' or 'C' key:  toggle screen clearing
        case 'C':					// to demonstrate 'trails'.
            if (isClear == 0) isClear = 1;
            else isClear = 0;
            break;
        case 't':			// REDUCE drag;  make velocity scale factor rise towards 1.0
            g_drag *= 1.0 / 0.995;
            if (g_drag > 1.0) g_drag = 1.0;	// don't allow drag to ADD energy!
            break;
        case 'T':			// INCREASE drag: make velocity scale factor a smaller fraction
            g_drag *= 0.995;
            break;
        case 'g':			// REDUCE gravity
            g_grav *= 0.99;		// shrink 1%
            break;
        case 'G':
            g_grav *= 1.0 / 0.98;	// grow by 2%
            break;
        case 'm':
            g_mass *= 0.98;		// reduce mass by 2%
            break;
        case 'M':
            g_mass *= 1.0 / 0.98;	// increase mass by 2%
            break;
        case 'R':  // HARD reset: position AND velocity.
            myRunMode = 0;			// RESET!
            s0[PART_XPOS] = 0.0;
            s0[PART_YPOS] = 0.0;
            s0[PART_ZPOS] = 0.0;
            s0[PART_XVEL] = INIT_VEL;
            s0[PART_YVEL] = INIT_VEL;
            s0[PART_ZVEL] = INIT_VEL;
            break;
        case 'r':		// 'SOFT' reset: boost velocity only.
            // don't change myRunMode
            if (g_show1 == 1) {
                part1Box.reset();
                break;
            }
            if (g_show3 == 1) {
                fluidParticleBox.reset(partSys1);
            }
        case 'z':
            if (g_show3 == 1) {
                fluidParticleBox.resetWind(partSys2);
            }
        case 'v':
        case 'V':
            // switch to a different solver:
            if (g_solver == 0) g_solver = 1;
            else g_solver = 0;
            break;
        case 'p':
        case 'P':			// toggle pause/run:
            if (myRunMode == 3) myRunMode = 1;		// if running, pause
            else myRunMode = 3;		// if paused, run.
            break;
        case ' ':			// space-bar: single-step
            myRunMode = 2;
            break;
        default:
            console.log('myKeyPress(): Ignored key: ' + myChar);
            break;
    }
}

VBObox0.prototype.displayMe = function () {
//==============================================================================
// Print current state of the particle system on the webpage:
    var recip = 1000.0 / timeStep;			// to report fractional seconds
    var solvType;												// convert solver number to text:
    if (g_solver == 0) solvType = 'Explicit--(unstable!)<br>';
    else solvType = 'Implicit--(stable)<br>';
    var bounceType;											// convert bounce number to text
    if (g_bounce == 0) bounceType = 'Velocity Reverse(no rest)<br>';
    else bounceType = 'Impulsive (will rest)<br>';
    var xvLimit = s0[PART_XVEL];							// find absolute values of s0[PART_XVEL]
    if (s0[PART_XVEL] < 0.0) xvLimit = -s0[PART_XVEL];
    var yvLimit = s0[PART_YVEL];							// find absolute values of s0[PART_YVEL]
    if (s0[PART_YVEL] < 0.0) yvLimit = -s0[PART_YVEL];
    var zvLimit = s0[PART_ZVEL];							// find absolute values of s0[PART_YVEL]
    if (s0[PART_ZVEL] < 0.0) zvLimit = -s0[PART_ZVEL];


    document.getElementById('KeyResult').innerHTML =
        '<b>Solver = </b>' + solvType +
        '<b>Bounce = </b>' + bounceType +
        '<b>drag = </b>' + g_drag.toFixed(5) +
        ', <b>grav = </b>' + g_grav.toFixed(5) +
        ' m/s^2; <b>yVel = +/-</b> ' + yvLimit.toFixed(5) +
        ' m/s; <b>xVel = +/-</b> ' + xvLimit.toFixed(5) +
        ' m/s; <b>zVel = +/-</b> ' + zvLimit.toFixed(5) +
        ' m/s;<br><b>timeStep = </b> 1/' + recip.toFixed(3) + ' sec<br>' +
        ' <b>stepCount: </b>' + g_stepCount.toFixed(3);
}


function onPlusButton() {
//==============================================================================
    INIT_VEL *= 1.2;		// increase
    console.log('Initial velocity: ' + INIT_VEL);
}

function onMinusButton() {
//==============================================================================
    INIT_VEL /= 1.2;		// shrink
    console.log('Initial velocity: ' + INIT_VEL);
}

/*
/
VBObox0.prototype.reload = function() {
//=============================================================================
// Over-write current values in the GPU inside our already-created VBO: use 
// gl.bufferSubData() call to re-transfer some or all of our Float32Array 
// contents to our VBO without changing any GPU memory allocations.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}

VBObox0.prototype.empty = function() {
//=============================================================================
// Remove/release all GPU resources used by this VBObox object, including any 
// shader programs, attributes, uniforms, textures, samplers or other claims on 
// GPU memory.  However, make sure this step is reversible by a call to 
// 'restoreMe()': be sure to retain all our Float32Array data, all values for 
// uniforms, all stride and offset values, etc.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}

VBObox0.prototype.restore = function() {
//=============================================================================
// Replace/restore all GPU resources used by this VBObox object, including any 
// shader programs, attributes, uniforms, textures, samplers or other claims on 
// GPU memory.  Use our retained Float32Array data, all values for  uniforms, 
// all stride and offset values, etc.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}
*/

//=============================================================================
//=============================================================================
function VBObox1() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox' object  that holds all data and fcns 
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate set of shaders.
    this.VERT_SRC =	//--------------------- VERTEX SHADER source code
        'precision mediump float;\n' +				// req'd in OpenGL ES if we use 'float'
        'uniform   int u_runMode1; \n' +					// particle system state:
        'uniform   mat4 u_MvpMat1; \n' +    //camera position
        'attribute vec3 a_Pos1;\n' +
        'attribute vec3 a_Colr1;\n' +
        'varying   vec4 v_Colr1; \n' +
        'attribute float a_diam1; \n' +
        'void main() {\n' +
        '  gl_PointSize = a_diam1;\n' +
        '	 gl_Position = u_MvpMat1 * (vec4(a_Pos1, 1)); \n' +// + u_ballShift1)

        // Let u_runMode determine particle color:
        'v_Colr1 = vec4(a_Colr1, 1.0); \n' +	// green: >3==run
        //	'		 } \n' +
        '} \n';

    this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code
        'precision mediump float;\n' +
        'varying vec4 v_Colr1; \n' +
        'void main() {\n' +
        '  float dist = distance(gl_PointCoord, vec2(0.5, 0.5)); \n' +
        '  	gl_FragColor = vec4((1.0-2.0*dist)*v_Colr1.rgb, 1.0);\n' +
        '}\n';


    this.vboLoc;										// Vertex Buffer Object location# on the GPU

    // bytes req'd for 1 vboContents array element;
    // (why? used to compute stride and offset
    // in bytes for vertexAttribPointer() calls)
    this.shaderLoc;									// Shader-program location # on the GPU, made
    // by compile/link of VERT_SRC and FRAG_SRC.
    //-------------------- Attribute locations in our shaders
    this.a_PosLoc;									// GPU location for 'a_Position' attribute
    this.a_ColrLoc;									// GPU location for 'a_Color' attribute

    //-------------------- Uniform locations &values in our shaders
    this.u_MvpMatLoc;								// GPU location for u_MvpMat uniform

    this.gndStart;


};


VBObox1.prototype.init = function (myGL) {
//=============================================================================
// Create, compile, link this VBObox object's shaders to an executable 'program'
// ready for use in the GPU.  Create and fill a Float32Array that holds all VBO 
// vertices' values; create a new VBO on the GPU and fill it with those values. 
// Find the GPU location of	all our shaders' attribute- and uniform-variables; 
// assign the correct portions of VBO contents as the data source for each 
// attribute, and transfer current values to the GPU for each uniform variable.
// (usually called only once, within main())
    this.s0 = new Float32Array(partCount * PART_MAXVAR);
    this.s1 = new Float32Array(partCount * PART_MAXVAR);
    this.FSIZE = this.s0.BYTES_PER_ELEMENT;	// memory needed to store an s0 array element.

    this.PartSys_init();
// Compile,link,upload shaders-------------------------------------------------
    this.shaderLoc = createProgram(myGL, this.VERT_SRC, this.FRAG_SRC);
    if (!this.shaderLoc) {
        console.log(this.constructor.name +
            '.init() failed to create executable Shaders on the GPU. Bye!');
        return;
    }
    // CUTE TRICK: we can print the NAME of this VBO object: tells us which one!
//  else{console.log('You called: '+ this.constructor.name + '.init() fcn!');}
    myGL.program = this.shaderLoc;		// (to match cuon-utils.js -- initShaders())

// Create VBO on GPU, fill it--------------------------------------------------
    this.vboLoc = myGL.createBuffer();
    if (!this.vboLoc) {
        console.log(this.constructor.name +
            '.init() failed to create VBO in GPU. Bye!');
        return;
    }
    // Specify the purpose of our newly-created VBO.  Your choices are:
    //	== "gl.ARRAY_BUFFER" : the VBO holds vertices, each made of attributes
    // (positions, colors, normals, etc), or
    //	== "gl.ELEMENT_ARRAY_BUFFER" : the VBO holds indices only; integer values
    // that each select one vertex from a vertex array stored in another VBO.
    myGL.bindBuffer(myGL.ARRAY_BUFFER,	// GLenum 'target' for this GPU buffer
        this.vboLoc);				// the ID# the GPU uses for this buffer.

    // Transfer data from JavaScript Float32Array object to the just-bound VBO.
    //  (Recall gl.bufferData() changes GPU's memory allocation: use
    //		gl.bufferSubData() to modify buffer contents without changing its size)
    //	The 'hint' helps GPU allocate its shared memory for best speed & efficiency
    //	(see OpenGL ES specification for more info).  Your choices are:
    //		--STATIC_DRAW is for vertex buffers rendered many times, but whose
    //				contents rarely or never change.
    //		--DYNAMIC_DRAW is for vertex buffers rendered many times, but whose
    //				contents may change often as our program runs.
    //		--STREAM_DRAW is for vertex buffers that are rendered a small number of
    // 			times and then discarded; for rapidly supplied & consumed VBOs.
    myGL.bufferData(gl.ARRAY_BUFFER, 			// GLenum target(same as 'bindBuffer()')
        this.s0, 		// JavaScript Float32Array
        gl.DYNAMIC_DRAW);			// Usage hint.

// Find & Set All Attributes:------------------------------
    // a) Get the GPU location for each attribute var used in our shaders:
    this.a_PosLoc = gl.getAttribLocation(this.shaderLoc, 'a_Pos1');
    if (this.a_PosLoc < 0) {
        console.log(this.constructor.name +
            '.init() Failed to get GPU location of attribute a_Pos1');
        return -1;	// error exit.
    }

    this.a_ColrLoc = myGL.getAttribLocation(this.shaderLoc, 'a_Colr1');
    if (this.a_ColrLoc < 0) {
        console.log(this.constructor.name +
            '.init() failed to get the GPU location of attribute a_Colr1');
        return -1;	// error exit.
    }
    this.a_diamID = myGL.getAttribLocation(this.shaderLoc, 'a_diam1');
    if (this.a_diamID < 0) {
        console.log(this.constructor.name +
            '.init() failed to get the GPU location of attribute a_diam1');
        return -1;	// error exit.
    }
    // b) Next, set up GPU to fill these attribute vars in our shader with
    // values pulled from the currently-bound VBO (see 'gl.bindBuffer()).
    // 	Here's how to use the almost-identical OpenGL version of this function:
    //		http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribPointer.xml )
    myGL.vertexAttribPointer(
        this.a_PosLoc,//index == ID# for the attribute var in your GLSL shaders;
        3,						// size == how many dimensions for this attribute: 1,2,3 or 4?
        gl.FLOAT,			// type == what data type did we use for those numbers?
        false,				// isNormalized == are these fixed-point values that we need
        //									normalize before use? true or false
        PART_MAXVAR * this.FSIZE,	// Stride == #bytes we must skip in the VBO to move from one
        //0,							// of our stored attributes to the next.  This is usually the
        // number of bytes used to store one complete vertex.  If set
        // to zero, the GPU gets attribute values sequentially from
        // VBO, starting at 'Offset'.
        // (Our vertex size in bytes: 4 floats for pos + 3 for color)
        0);						// Offset == how many bytes from START of buffer to the first
    // value we will actually use?  (We start with position).
    gl.vertexAttribPointer(this.a_ColrLoc, 3, gl.FLOAT, false,
        PART_MAXVAR * this.FSIZE, PART_R * this.FSIZE);
    myGL.vertexAttribPointer(this.a_diamID, 1, gl.FLOAT, false,
        PART_MAXVAR * this.FSIZE, PART_DIAM * this.FSIZE);
    // c) Enable this assignment of the attribute to its' VBO source:
    myGL.enableVertexAttribArray(this.a_PosLoc);
    myGL.enableVertexAttribArray(this.a_ColrLoc);
    myGL.enableVertexAttribArray(this.a_diamID);
// Find All Uniforms:--------------------------------
//Get GPU storage location for each uniform var used in our shader programs: 
    this.u_MvpMatLoc = myGL.getUniformLocation(this.shaderLoc, 'u_MvpMat1');
    if (!this.u_MvpMatLoc) {
        console.log(this.constructor.name +
            '.init() failed to get GPU location for u_MvpMat1 uniform');
        return;
    }

    this.u_runModeID = myGL.getUniformLocation(this.shaderLoc, 'u_runMode1');
    if (!this.u_runModeID) {
        console.log('Failed to get u_runModeID variable location');
        return;
    }

}

VBObox1.prototype.PartSys_init = function () {
//==============================================================================
// set initial values of all particle-system state.
// sel==0 for 'hard' reset (user pressed 'R') that inits entire state-vector
// sel==1 to 'soft' reset (user pressed 'r') that only adds velocity to all 
//						particles

    var doit = 1;
    for (var i = 0; i < partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        var xcyc = roundRand3D();
        if (doit == 1) {
            console.log('xc,yc= ' + xcyc[0] + ', ' + xcyc[1]);
            doit = 0;
        }
        this.s0[pOff + PART_XPOS] = 0.2 + 0.2 * xcyc[0];		// 0.0 <= randomRound() < 1.0
        this.s0[pOff + PART_YPOS] = 0.2 + 0.2 * xcyc[1];
        this.s0[pOff + PART_ZPOS] = 0.2 + 0.2 * xcyc[2];
        xcyc = roundRand3D();
        this.s0[pOff + PART_XVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[0]);
        this.s0[pOff + PART_YVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[1]);
        this.s0[pOff + PART_ZVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[2]);
        this.s0[pOff + PART_X_FTOT] = 0.0;
        this.s0[pOff + PART_Y_FTOT] = 0.0;
        this.s0[pOff + PART_Z_FTOT] = 0.0;
        this.s0[pOff + PART_R] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_G] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_B] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_MASS] = 0.9 + 0.2 * Math.random();
        this.s0[pOff + PART_DIAM] = 1.0 + 10.0 * Math.random();
        this.s0[pOff + PART_RENDMODE] = Math.floor(4.0 * Math.random()); // 0,1,2 or 3.
    }
}
VBObox1.prototype.adjust = function (myGL) {
//=============================================================================
// Update the GPU to newer, current values we now store for 'uniform' vars on 
// the GPU; and (if needed) update each attribute's stride and offset in VBO.

    myGL.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    // -and-linked executable shader program.
    for (i = 0; i < partCount; i++)
        if (myRunMode > 1) {											// 0=reset; 1= pause; 2=step; 3=run
            if (myRunMode == 2) myRunMode = 1;				// (if 2, do just one step and pause.)
            //=YES!=========================================
            // Make our 'bouncy-ball' move forward by one timestep, but now the 's' key
            // will select which kind of solver to use:
            if (g_solver == 0) {
                //-----------------------------------------------------------------------
                // EXPLICIT or 'forward time' solver, as found in bouncyBall03.01BAD and
                // bouncyBall04.01badMKS.  CAREFUL! this solver adds energy -- not stable
                // for many particle system settings!
                // This solver looks quite sensible and logical.  Formally, it's an
                //	explicit or 'forward-time' solver known as the Euler method:
                //			Use the current velocity ('s0dot') to move forward by
                //			one timestep: s1 = s0 + s0dot*h, and
                //		-- Compute the new velocity (e.g. s1dot) too: apply gravity & drag.
                //		-- Then apply constraints: check to see if new position (s1)
                //			is outside our floor, ceiling, or walls, and if new velocity
                //			will move us further in the wrong direction. If so, reverse it!
                // CAREFUL! must convert timeStep from milliseconds to seconds!
                xposPrev = this.s0[i * PART_MAXVAR + PART_XPOS];			// SAVE these values before we update them.
                xvelPrev = this.s0[i * PART_MAXVAR + PART_XVEL];			// (for use in constraint-applying code below).
                yposPrev = this.s0[i * PART_MAXVAR + PART_YPOS];
                yvelPrev = this.s0[i * PART_MAXVAR + PART_YVEL];
                zposPrev = this.s0[i * PART_MAXVAR + PART_ZPOS];
                zvelPrev = this.s0[i * PART_MAXVAR + PART_ZVEL];

                //------------------
                // Compute new position from current position, current velocity, & timestep
                this.s0[i * PART_MAXVAR + PART_XPOS] += this.s0[i * PART_MAXVAR + PART_XVEL] * (timeStep * 0.001);
                this.s0[i * PART_MAXVAR + PART_YPOS] += this.s0[i * PART_MAXVAR + PART_YVEL] * (timeStep * 0.001);
                this.s0[i * PART_MAXVAR + PART_ZPOS] += this.s0[i * PART_MAXVAR + PART_ZVEL] * (timeStep * 0.001);
                // -- apply acceleration due to gravity to current velocity:
                // 					 this.s0[PART_ZVEL] -= (accel. due to gravity)*(timestep in seconds)
                //									 -= (9.832 meters/sec^2) * (timeStep/1000.0);
                this.s0[i * PART_MAXVAR + PART_ZVEL] -= g_grav * (timeStep * 0.001);
                // -- apply drag: attenuate current velocity:
                this.s0[i * PART_MAXVAR + PART_XVEL] *= g_drag;
                this.s0[i * PART_MAXVAR + PART_YVEL] *= g_drag;
                this.s0[i * PART_MAXVAR + PART_ZVEL] *= g_drag;

            } else if (g_solver == 1) {
                //------------------------------------------------------------------------
                // IMPLICIT or 'reverse time' solver, as found in bouncyBall04.goodMKS;
                // This category of solver is often better, more stable, but lossy.
                // -- apply acceleration due to gravity to current velocity:
                //				  s0[PART_YVEL] -= (accel. due to gravity)*(timestep in seconds)
                //                  -= (9.832 meters/sec^2) * (timeStep/1000.0);
                xposPrev = this.s0[i * PART_MAXVAR + PART_XPOS];			// SAVE these values before we update them.
                xvelPrev = this.s0[i * PART_MAXVAR + PART_XVEL];			// (for use in constraint-applying code below).
                yposPrev = this.s0[i * PART_MAXVAR + PART_YPOS];
                yvelPrev = this.s0[i * PART_MAXVAR + PART_YVEL];
                zposPrev = this.s0[i * PART_MAXVAR + PART_ZPOS];
                zvelPrev = this.s0[i * PART_MAXVAR + PART_ZVEL];
                //-------------------
                this.s0[i * PART_MAXVAR + PART_ZVEL] -= g_grav * (timeStep * 0.001);
                // -- apply drag: attenuate current velocity:
                this.s0[i * PART_MAXVAR + PART_XVEL] *= g_drag;
                this.s0[i * PART_MAXVAR + PART_YVEL] *= g_drag;
                this.s0[i * PART_MAXVAR + PART_ZVEL] *= g_drag;
                // -- move our particle using current velocity:
                // CAREFUL! must convert timeStep from milliseconds to seconds!
                this.s0[i * PART_MAXVAR + PART_XPOS] += this.s0[i * PART_MAXVAR + PART_XVEL] * (timeStep * 0.001);
                this.s0[i * PART_MAXVAR + PART_YPOS] += this.s0[i * PART_MAXVAR + PART_YVEL] * (timeStep * 0.001);
                this.s0[i * PART_MAXVAR + PART_ZPOS] += this.s0[i * PART_MAXVAR + PART_ZVEL] * (timeStep * 0.001);
                // What's the result of this rearrangement?
                //	IT WORKS BEAUTIFULLY! much more stable much more often...
            } else {
                console.log('?!?! unknown solver: g_solver==' + g_solver);
                return;
            }

            //==========================================================================
            // CONSTRAINTS -- 'bounce' our ball off floor & walls at (0,0), (2.0, 2.0):
            // where g_bounce selects constraint type:
            // ==0 for simple velocity-reversal, as in all previous versions
            // ==1 for Chapter 7's collision resolution method, which uses an 'impulse'
            //      to cancel any velocity boost caused by falling below the floor.
            if (g_bounce == 0) { //--------------------------------------------------------
                if (this.s0[i * PART_MAXVAR + PART_XPOS] < 0.0 && this.s0[i * PART_MAXVAR + PART_XVEL] < 0.0			// simple velocity-reversal
                ) {		// bounce on left wall.
                    this.s0[i * PART_MAXVAR + PART_XVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_XVEL];
                } else if (this.s0[i * PART_MAXVAR + PART_XPOS] > 2.0 && this.s0[i * PART_MAXVAR + PART_XVEL] > 0.0
                ) {		// bounce on right wall
                    this.s0[i * PART_MAXVAR + PART_XVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_XVEL];
                }
                if (this.s0[i * PART_MAXVAR + PART_YPOS] < 0.0 && this.s0[i * PART_MAXVAR + PART_YVEL] < 0.0
                ) {		// bounce on floor
                    this.s0[i * PART_MAXVAR + PART_YVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_YVEL];
                } else if (this.s0[i * PART_MAXVAR + PART_YPOS] > 2.0 && this.s0[i * PART_MAXVAR + PART_YVEL] > 0.0
                ) {		// bounce on ceiling
                    this.s0[i * PART_MAXVAR + PART_YVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_YVEL];
                }
                if (this.s0[i * PART_MAXVAR + PART_ZPOS] < 0 && this.s0[i * PART_MAXVAR + PART_ZVEL] < 0.0			// simple velocity-reversal
                ) {		// bounce on left wall.
                    this.s0[i * PART_MAXVAR + PART_ZVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_ZVEL];
                } else if (this.s0[i * PART_MAXVAR + PART_ZPOS] > 2.0 && this.s0[i * PART_MAXVAR + PART_ZVEL] > 0.0
                ) {		// bounce on right wall
                    this.s0[i * PART_MAXVAR + PART_ZVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_ZVEL];
                }
                //  -- hard limit on 'floor' keeps y position >= 0;
                if (this.s0[i * PART_MAXVAR + PART_ZPOS] < 0.0) this.s0[i * PART_MAXVAR + PART_ZPOS] = 0.0;
            } else if (g_bounce == 1) { //---------------------------------------------------------------------------
                if (this.s0[i * PART_MAXVAR + PART_XPOS] < 0.0 && this.s0[i * PART_MAXVAR + PART_XVEL] < 0.0 // collision!  left wall...
                ) {		// bounce on left wall.

                    this.s0[i * PART_MAXVAR + PART_XPOS] = 0.0;					// 1) resolve contact: put particle at wall.

                    this.s0[i * PART_MAXVAR + PART_XVEL] = xvelPrev;			// we had a the START of the timestep.
                    // (NOTE: statistically, hitting the wall is equally probable at any
                    // time during the timestep, so the 'expected value' of collision is at
                    // the timestep's midpoint.  THUS removing HALF the new velocity during
                    // the timestep would create errors with a statistical mean of zero.
                    //
                    // 		Unwittingly, we have already created that result!
                    //============================================================
                    // For simplicity, assume our timestep's erroneous velocity change
                    // was the result of constant acceleration (e.g. result of constant
                    // gravity acting constant mass, plus constant drag force, etc).  If the
                    // ball 'bounces' (reverses velocity) exactly halfway through the
                    // timestep, then at the statistical 'expected value' for collision
                    // time, then the constant force that acts to INCREASE velocity in one
                    // half-timestep will act to DECREASE velocity in the other half
                    // timestep by exactly the same amount -- and thus removes ALL the
                    // velocity added by constant force during the timestep.)
                    this.s0[i * PART_MAXVAR + PART_XVEL] *= g_drag;			// **BUT** velocity during our timestep is STILL
                    // reduced by drag (and any other forces
                    // proportional to velocity, and thus not
                    // cancelled by 'bounce' at timestep's midpoint)
                    // 3) BOUNCE:
                    //reversed velocity*coeff-of-restitution.
                    // ATTENTION! VERY SUBTLE PROBLEM HERE! ------------------------------
                    //Balls with tiny, near-zero velocities (e.g. ball nearly at rest on
                    // floor) can easily reverse sign between 'previous' and 'now'
                    // timesteps, even for negligible forces.  Put another way:
                    // Step 2), our 'repair' attempt that removes all erroneous x velocity,
                    // has CHANGED the 'now' ball velocity, and MAY have changed its sign as
                    // well,  especially when the ball is nearly at rest. SUBTLE: THUS we
                    // need a velocity-sign test here that ensures the 'bounce' step will
                    // always send the ball outwards, away from its wall or floor collision.
                    if (this.s0[i * PART_MAXVAR + PART_XVEL] < 0.0) this.s0[i * PART_MAXVAR + PART_XVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_XVEL]; // no sign change--bounce!
                    else this.s0[i * PART_MAXVAR + PART_XVEL] = g_resti * this.s0[i * PART_MAXVAR + PART_XVEL];			// sign changed-- don't need another.

                } else if (this.s0[i * PART_MAXVAR + PART_XPOS] > 2.0 && this.s0[i * PART_MAXVAR + PART_XVEL] > 0.0		// collision! right wall...
                ) {		// bounce on right wall

                    this.s0[i * PART_MAXVAR + PART_XPOS] = 2.0;
                    this.s0[i * PART_MAXVAR + PART_XVEL] = xvelPrev;			// velocity we had at the start of timestep;
                    this.s0[i * PART_MAXVAR + PART_XVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    // 	that still apply during this timestep).

                    if (this.s0[i * PART_MAXVAR + PART_XVEL] > 0.0) this.s0[i * PART_MAXVAR + PART_XVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_XVEL]; // no sign change--bounce!
                    else this.s0[i * PART_MAXVAR + PART_XVEL] = g_resti * this.s0[i * PART_MAXVAR + PART_XVEL];			// sign changed-- don't need another.

                }
                if (this.s0[i * PART_MAXVAR + PART_YPOS] < 0.0 && this.s0[i * PART_MAXVAR + PART_YVEL] < 0.0		// collision! left wall...
                ) {		// bounce on floor

                    this.s0[i * PART_MAXVAR + PART_YPOS] = 0.0;					// 1) resolve contact: put particle at wall.
                    // 2) remove all y velocity gained from forces as
                    // ball moved thru floor in this timestep. HOW?
                    // Assume ball reached floor at START of
                    // the timestep, thus: return to the orig.
                    this.s0[i * PART_MAXVAR + PART_YVEL] = yvelPrev;			// velocity we had at the start of timestep;
                    this.s0[i * PART_MAXVAR + PART_YVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    // 	that still apply during this timestep).

                    if (this.s0[i * PART_MAXVAR + PART_YVEL] < 0.0) this.s0[i * PART_MAXVAR + PART_YVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_YVEL]; // no sign change--bounce!
                    else this.s0[i * PART_MAXVAR + PART_YVEL] = g_resti * this.s0[i * PART_MAXVAR + PART_YVEL];			// sign changed-- don't need another.

                } else if (this.s0[i * PART_MAXVAR + PART_YPOS] > 2.0 && this.s0[i * PART_MAXVAR + PART_YVEL] > 0.0 		// collision! front wall...
                ) {		// bounce on ceiling

                    this.s0[i * PART_MAXVAR + PART_YPOS] = 2.0;					// 1) resolve contact: put particle at wall.
                    // 2) remove all y velocity gained from forces as
                    // ball moved thru ceiling in this timestep. HOW?
                    // Assume ball reached ceiling at START of
                    // the timestep, thus: return to the orig.
                    this.s0[i * PART_MAXVAR + PART_YVEL] = yvelPrev;			// velocity we had at the start of timestep;
                    this.s0[i * PART_MAXVAR + PART_YVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    // 	that still apply during this timestep),

                    if (this.s0[i * PART_MAXVAR + PART_YVEL] > 0.0) this.s0[i * PART_MAXVAR + PART_YVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_YVEL]; // no sign change--bounce!
                    else this.s0[i * PART_MAXVAR + PART_YVEL] = g_resti * this.s0[i * PART_MAXVAR + PART_YVEL];			// sign changed-- don't need another.

                }
                //	console.log('z = ' + s0[PART_ZPOS] + '  zVel = ' +s0[PART_ZVEL]);
                if (this.s0[i * PART_MAXVAR + PART_ZPOS] < 0.0 && this.s0[i * PART_MAXVAR + PART_ZVEL] < 0.0		// collision! left wall...
                ) {		// bounce on floor

                    this.s0[i * PART_MAXVAR + PART_ZPOS] = 0.0;					// 1) resolve contact: put particle at wall.
                    // 2) remove all y velocity gained from forces as
                    // ball moved thru floor in this timestep. HOW?
                    // Assume ball reached floor at START of
                    // the timestep, thus: return to the orig.
                    this.s0[i * PART_MAXVAR + PART_ZVEL] = zvelPrev;			// velocity we had at the start of timestep;
                    this.s0[i * PART_MAXVAR + PART_ZVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    // 	that still apply during this timestep).

                    if (this.s0[i * PART_MAXVAR + PART_ZVEL] < 0.0) this.s0[i * PART_MAXVAR + PART_ZVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_ZVEL]; // no sign change--bounce!
                    else this.s0[i * PART_MAXVAR + PART_ZVEL] = g_resti * this.s0[i * PART_MAXVAR + PART_ZVEL];			// sign changed-- don't need another.

                } else if (this.s0[i * PART_MAXVAR + PART_ZPOS] > 2.0 && this.s0[i * PART_MAXVAR + PART_ZVEL] > 0.0 		// collision! front wall...
                ) {		// bounce on ceiling

                    this.s0[i * PART_MAXVAR + PART_ZPOS] = 2.0;					// 1) resolve contact: put particle at wall.
                    // 2) remove all y velocity gained from forces as
                    // ball moved thru ceiling in this timestep. HOW?
                    // Assume ball reached ceiling at START of
                    // the timestep, thus: return to the orig.
                    this.s0[i * PART_MAXVAR + PART_ZVEL] = zvelPrev;			// velocity we had at the start of timestep;
                    this.s0[i * PART_MAXVAR + PART_ZVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    // 	that still apply during this timestep),

                    if (this.s0[i * PART_MAXVAR + PART_ZVEL] > 0.0) this.s0[i * PART_MAXVAR + PART_ZVEL] = -g_resti * this.s0[i * PART_MAXVAR + PART_ZVEL]; // no sign change--bounce!
                    else this.s0[i * PART_MAXVAR + PART_ZVEL] = g_resti * this.s0[i * PART_MAXVAR + PART_ZVEL];			// sign changed-- don't need another.

                }
            } else {
                console.log('?!?! unknown constraint: g_bounce==' + g_bounce);
                return;
            }
            //============================================
        }

    // Adjust values for our uniforms,

    var vpAspect = ((myGL.drawingBufferWidth) /			// On-screen aspect ratio for
        myGL.drawingBufferHeight);		// this camera: width/height.

    MvpMat.setIdentity();	// rotate drawing axes,
    MvpMat.perspective(45, vpAspect, 1, +200);							// then translate them.
    MvpMat.lookAt(xCamPos, yCamPos, zCamPos, xCamPos + Math.cos(Math.PI * (currentAngle / 180)),
        yCamPos + Math.sin(Math.PI * (currentAngle / 180)), zCamPos + zlook, 0.0, 0, 1);   // UP vector.

    //  Transfer new uniforms' values to the GPU:-------------
    // Send  new 'MvpMat' values to the GPU's 'u_MvpMat1' uniform:
    myGL.uniformMatrix4fv(this.u_MvpMatLoc,	// GPU location of the uniform
        false, 				// use matrix transpose instead?
        MvpMat.elements);	// send data from Javascript.
    myGL.uniform1i(this.u_runModeID, myRunMode);		// run/step/pause the particle system

}

VBObox1.prototype.PartSys_render = function (myGL) {
//=============================================================================
// Send commands to GPU to select and render current VBObox contents.  	
    myGL.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    // -and-linked executable shader program.
//------CAREFUL! RE-BIND YOUR VBO AND RE-ASSIGN SHADER ATTRIBUTES!-------------
//		Each call to useProgram() reconfigures the GPU's processors & data paths 
// for efficient SIMD execution of the newly-selected shader program. While the 
// 'old' shader program's attributes and uniforms remain at their same memory 
// locations, starting the new shader program invalidates the old data paths 
// that connected these attributes to the VBOs in memory that supplied their 
// values. When we call useProgram() to return to our 'old' shader program, we 
// must re-establish those data-paths between shader attributes and VBOs, even 
// if those attributes, VBOs, and locations have not changed!
//		Thus after each useProgram() call, we must:
// a)--call bindBuffer() again to re-bind each VBO that our shader will use, &
// b)--call vertexAttribPointer() again for each attribute in our new shader
//		program, to re-connect the data-path(s) from bound VBO(s) to attribute(s):
// c)--call enableVertexAttribArray() to enable use of those data paths.
//----------------------------------------------------
    // a) Re-set the GPU's currently 'bound' vbo buffer;
    myGL.bindBuffer(myGL.ARRAY_BUFFER,	// GLenum 'target' for this GPU buffer
        this.vboLoc);			// the ID# the GPU uses for this buffer.
    // (Here's how to use the almost-identical OpenGL version of this function:
    //		http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribPointer.xml )
    // b) Re-connect data paths from VBO to each shader attribute:
    myGL.vertexAttribPointer(this.a_PosLoc, 3, myGL.FLOAT, false,
        PART_MAXVAR * this.FSIZE, 0);		// stride, offset
    myGL.vertexAttribPointer(this.a_ColrLoc, 3, myGL.FLOAT, false,
        PART_MAXVAR * this.FSIZE, PART_R * this.FSIZE); // stride, offset
    myGL.vertexAttribPointer(this.a_diamID, 1, myGL.FLOAT, false,
        PART_MAXVAR * this.FSIZE, PART_DIAM * this.FSIZE); // stride, offset
    // c) enable the newly-re-assigned attributes:
    myGL.enableVertexAttribArray(this.a_PosLoc);
    myGL.enableVertexAttribArray(this.a_ColrLoc);
    myGL.enableVertexAttribArray(this.a_diamID);

    // ----------------------------Draw the contents of the currently-bound VBO:


    myGL.bufferSubData(myGL.ARRAY_BUFFER, 0, this.s0);
    myGL.drawArrays(myGL.POINTS, 0, partCount);
}
/*
/
VBObox1.prototype.reload = function() {
//=============================================================================
// Over-write current values in the GPU for our already-created VBO: use 
// gl.bufferSubData() call to re-transfer some or all of our Float32Array 
// contents to our VBO without changing any GPU memory allocations.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}

VBObox1.prototype.empty = function() {
//=============================================================================
// Remove/release all GPU resources used by this VBObox object, including any 
// shader programs, attributes, uniforms, textures, samplers or other claims on 
// GPU memory.  However, make sure this step is reversible by a call to 
// 'restoreMe()': be sure to retain all our Float32Array data, all values for 
// uniforms, all stride and offset values, etc.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}

VBObox1.prototype.restore = function() {
//=============================================================================
// Replace/restore all GPU resources used by this VBObox object, including any 
// shader programs, attributes, uniforms, textures, samplers or other claims on 
// GPU memory.  Use our retained Float32Array data, all values for  uniforms, 
// all stride and offset values, etc.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}
*/

//=============================================================================
//=============================================================================
function VBObox2() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox' object  that holds all data and fcns 
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate set of shaders.

    this.VERT_SRC =	//--------------------- VERTEX SHADER source code
        'precision highp float;\n' +				// req'd in OpenGL ES if we use 'float'
        //
        'uniform mat4 u_ModelMatrix;\n' +
        'attribute vec4 a_Position;\n' +
        'attribute vec3 a_Color;\n' +
        'attribute float a_PtSize; \n' +
        'varying vec3 v_Colr2;\n' +
        //
        'void main() {\n' +
        '  gl_PointSize = a_PtSize;\n' +
        '  gl_Position = u_ModelMatrix * a_Position;\n' +
        '	 v_Colr2 = a_Color;\n' +
        ' }\n';

    this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code
        'precision mediump float;\n' +
        'varying vec3 v_Colr2;\n' +
        'void main() {\n' +
        '  gl_FragColor = vec4(v_Colr2, 1.0);\n' +
        '}\n';

    this.vboContents = //---------------------
        new Float32Array([					// Array of vertex attribute values we will
            // transfer to GPU's vertex buffer object (VBO)
            // 1 vertex per line: pos x,y,z,w;   color; r,g,b;   point-size;
            -0.3, 0.7, 0.0, 1.0, 1.0, 0.3, 0.3, 13.0,
            -0.3, -0.3, 0.0, 1.0, 0.3, 1.0, 0.3, 20.0,
            0.3, 0, 0.0, 1.0, 0.3, 0.3, 1.0, 7.0,
        ]);
    this.vboVerts = 3;							// # of vertices held in 'vboContents' array;
    this.vboLoc;										// Vertex Buffer Object location# on the GPU
    this.FSIZE = this.vboContents.BYTES_PER_ELEMENT;
    // bytes req'd for 1 array element;
    // (why? used to compute stride and offset
    // in bytes for vertexAttribPointer() calls)
    this.shaderLoc;									// Shader-program location # on the GPU, made
    // by compile/link of VERT_SRC and FRAG_SRC.
    //-------------------- Attribute locations in our shaders
    this.a_PositionLoc;							// GPU location: shader 'a_Position' attribute
    this.a_ColorLoc;								// GPU location: shader 'a_Color' attribute
    this.a_PtSizeLoc;								// GPU location: shader 'a_PtSize' attribute
    //-------------------- Uniform locations &values in our shaders
    this.ModelMatrix = new Matrix4();	// Transforms CVV axes to model axes.
    this.u_ModelMatrixLoc;						// GPU location for u_ModelMat uniform
};


VBObox2.prototype.init = function (myGL) {
//=============================================================================
// Create, compile, link this VBObox object's shaders to an executable 'program'
// ready for use in the GPU.  Create and fill a Float32Array that holds all VBO 
// vertices' values; create a new VBO on the GPU and fill it with those values. 
// Find the GPU location of	all our shaders' attribute- and uniform-variables; 
// assign the correct portions of VBO contents as the data source for each 
// attribute, and transfer current values to the GPU for each uniform variable.
// (usually called only once, within main()) 

    this.shaderLoc = createProgram(myGL, this.VERT_SRC, this.FRAG_SRC);
    if (!this.shaderLoc) {
        console.log(this.constructor.name +
            '.init() failed to create executable Shaders on the GPU. Bye!');
        return;
    }
    this.vboLoc = myGL.createBuffer();
    if (!this.vboLoc) {
        console.log(this.constructor.name +
            '.init() failed to create VBO in GPU. Bye!');
        return;
    }
    // Specify the purpose of our newly-created VBO.  Your choices are:
    //	== "gl.ARRAY_BUFFER" : the VBO holds vertices, each made of attributes
    // (positions, colors, normals, etc), or
    //	== "gl.ELEMENT_ARRAY_BUFFER" : the VBO holds indices only; integer values
    // that each select one vertex from a vertex array stored in another VBO.
    myGL.bindBuffer(myGL.ARRAY_BUFFER,	// GLenum 'target' for this GPU buffer
        this.vboLoc);				// the ID# the GPU uses for this buffer.

    // Transfer data from our JavaScript Float32Array object to the just-bound VBO.
    //  (Recall gl.bufferData() changes GPU's memory allocation: use
    //		gl.bufferSubData() to modify buffer contents without changing its size)
    //	The 'hint' helps GPU allocate its shared memory for best speed & efficiency
    //	(see OpenGL ES specification for more info).  Your choices are:
    //		--STATIC_DRAW is for vertex buffers rendered many times, but whose
    //				contents rarely or never change.
    //		--DYNAMIC_DRAW is for vertex buffers rendered many times, but whose
    //				contents may change often as our program runs.
    //		--STREAM_DRAW is for vertex buffers that are rendered a small number of
    // 			times and then discarded; for rapidly supplied & consumed VBOs.
    myGL.bufferData(gl.ARRAY_BUFFER, 			// GLenum target(same as 'bindBuffer()')
        this.vboContents, 		// JavaScript Float32Array
        gl.STATIC_DRAW);			// Usage hint.

// Find & Set All Attributes:------------------------------
    // a) Get the GPU location for each attribute var used in our shaders:
    this.a_PositionLoc = gl.getAttribLocation(this.shaderLoc, 'a_Position');
    if (this.a_PositionLoc < 0) {
        console.log(this.constructor.name +
            '.init() Failed to get GPU location of attribute a_Position');
        return -1;	// error exit.
    }
    this.a_ColorLoc = myGL.getAttribLocation(this.shaderLoc, 'a_Color');
    if (this.a_ColorLoc < 0) {
        console.log(this.constructor.name +
            '.init() failed to get the GPU location of attribute a_Color');
        return -1;	// error exit.
    }
    // NEW! a_PtSize' attribute values are stored only in VBO2, not VBO1:
    this.a_PtSizeLoc = gl.getAttribLocation(this.shaderLoc, 'a_PtSize');
    if (this.a_PtSizeLoc < 0) {
        console.log(this.constructor.name +
            '.init() failed to get the GPU location of attribute a_PtSize');
        return -1;	// error exit.
    }
    // b) Next, set up GPU to fill these attribute vars in our shader with
    // values pulled from the currently-bound VBO (see 'gl.bindBuffer()).
    // 	Here's how to use the almost-identical OpenGL version of this function:
    //		http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribPointer.xml )
    myGL.vertexAttribPointer(
        this.a_PositionLoc,//index == ID# for the attribute var in GLSL shader pgm;
        4,						// size == how many dimensions for this attribute: 1,2,3 or 4?
        myGL.FLOAT,		// type == what data type did we use for those numbers?
        false,				// isNormalized == are these fixed-point values that we need
        //									normalize before use? true or false
        8 * this.FSIZE,	// Stride == #bytes we must skip in the VBO to move from one
        // of our stored attributes to the next.  This is usually the
        // number of bytes used to store one complete vertex.  If set
        // to zero, the GPU gets attribute values sequentially from
        // VBO, starting at 'Offset'.
        // (Vertex size in bytes: 4 floats for pos; 3 color; 1 PtSize)
        0);						// Offset == how many bytes from START of buffer to the first
    // value we will actually use?  (We start with position).
    myGL.vertexAttribPointer(this.a_ColorLoc, 3, myGL.FLOAT, false,
        8 * this.FSIZE, 			// stride for VBO2 (different from VBO1!)
        4 * this.FSIZE);			// offset: skip the 1st 4 floats.
    myGL.vertexAttribPointer(this.a_PtSizeLoc, 1, myGL.FLOAT, false,
        8 * this.FSIZE,		// stride for VBO2 (different from VBO1!)
        7 * this.FSIZE);		// offset: skip the 1st 7 floats.

    // c) Enable this assignment of the attribute to its' VBO source:
    myGL.enableVertexAttribArray(this.a_PositionLoc);
    myGL.enableVertexAttribArray(this.a_ColorLoc);
    myGL.enableVertexAttribArray(this.a_PtSizeLoc);

// Find All Uniforms:--------------------------------
//Get GPU storage location for each uniform var used in our shader programs: 
    this.u_ModelMatrixLoc = gl.getUniformLocation(this.shaderLoc, 'u_ModelMatrix');
    if (!this.u_ModelMatrixLoc) {
        console.log(this.constructor.name +
            '.init() failed to get GPU location for u_ModelMatrix uniform');
        return;
    }
}

VBObox2.prototype.adjust = function (myGL) {
//=============================================================================
// Update the GPU to newer, current values we now store for 'uniform' vars on 
// the GPU; and (if needed) update each attribute's stride and offset in VBO.
    myGL.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    // -and-linked executable shader program.
    // Adjust values for our uniforms: -----------------------
    this.ModelMatrix.setRotate(-g_currentAngle * 2, 0, 0, 1);	// -spin drawing axes,
    this.ModelMatrix.translate(-0.35, -0.25, 0);							// then translate them.
    //  Transfer new uniforms' values to the GPU:-------------
    // Send  new 'ModelMat' values to the GPU's 'u_ModelMat1' uniform:
    myGL.uniformMatrix4fv(this.u_ModelMatrixLoc,	// GPU location of the uniform
        false, 										// use matrix transpose instead?
        this.ModelMatrix.elements);	// send data from Javascript.
}

VBObox2.prototype.draw = function (myGL) {
//=============================================================================
// Send commands to GPU to select and render current VBObox contents.

    myGL.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    // -and-linked executable shader program.
//------CAREFUL! RE-BIND YOUR VBO AND RE-ASSIGN SHADER ATTRIBUTES!-------------
//		Each call to useProgram() reconfigures the GPU's processors & data paths 
// for efficient SIMD execution of the newly-selected shader program. While the 
// 'old' shader program's attributes and uniforms remain at their same memory 
// locations, starting the new shader program invalidates the old data paths 
// that connected these attributes to the VBOs in memory that supplied their 
// values. When we call useProgram() to return to our 'old' shader program, we 
// must re-establish those data-paths between shader attributes and VBOs, even 
// if those attributes, VBOs, and locations have not changed!
//		Thus after each useProgram() call, we must:
// a)--call bindBuffer() again to re-bind each VBO that our shader will use, &
// b)--call vertexAttribPointer() again for each attribute in our new shader
//		program, to re-connect the data-path(s) from bound VBO(s) to attribute(s):
// c)--call enableVertexAttribArray() to enable use of those data paths.
//----------------------------------------------------
    // a) Re-set the GPU's currently 'bound' vbo buffer;
    myGL.bindBuffer(myGL.ARRAY_BUFFER,	// GLenum 'target' for this GPU buffer
        this.vboLoc);			// the ID# the GPU uses for this buffer.
    // (Here's how to use the almost-identical OpenGL version of this function:
    //		http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribPointer.xml )
    //b) Re-connect data paths from VBO to each shader attribute:
    myGL.vertexAttribPointer(this.a_PositionLoc, 4, myGL.FLOAT, false,
        8 * this.FSIZE, 0);							// stride, offset
    myGL.vertexAttribPointer(this.a_ColorLoc, 3, myGL.FLOAT, false,
        8 * this.FSIZE, 4 * this.FSIZE);	// stride, offset
    myGL.vertexAttribPointer(this.a_PtSizeLoc, 1, myGL.FLOAT, false,
        8 * this.FSIZE, 7 * this.FSIZE);	// stride, offset
    // c) Re-Enable use of the data path for each attribute:
    myGL.enableVertexAttribArray(this.a_PositionLoc);
    myGL.enableVertexAttribArray(this.a_ColorLoc);
    myGL.enableVertexAttribArray(this.a_PtSizeLoc);
    // ----------------------------Draw the contents of the currently-bound VBO:
    myGL.drawArrays(myGL.POINTS, 		// select the drawing primitive to draw,
        0, 								// location of 1st vertex to draw;
        this.vboVerts);		// number of vertices to draw on-screen.

    myGL.drawArrays(myGL.LINE_LOOP, // draw lines between verts too
        0,
        this.vboVerts);
}
/*
/
VBObox2.prototype.reload = function() {
//=============================================================================
// Over-write current values in the GPU for our already-created VBO: use 
// gl.bufferSubData() call to re-transfer some or all of our Float32Array 
// contents to our VBO without changing any GPU memory allocations.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}

VBObox2.prototype.empty = function() {
//=============================================================================
// Remove/release all GPU resources used by this VBObox object, including any 
// shader programs, attributes, uniforms, textures, samplers or other claims on 
// GPU memory.  However, make sure this step is reversible by a call to 
// 'restoreMe()': be sure to retain all our Float32Array data, all values for 
// uniforms, all stride and offset values, etc.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}

VBObox2.prototype.restore = function() {
//=============================================================================
// Replace/restore all GPU resources used by this VBObox object, including any 
// shader programs, attributes, uniforms, textures, samplers or other claims on 
// GPU memory.  Use our retained Float32Array data, all values for  uniforms, 
// all stride and offset values, etc.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}
*/

/*
//=============================================================================
//=============================================================================
//=============================================================================
function animate() {
//==============================================================================  
// How much time passed since we last updated the 'canvas' screen elements?
  var now = Date.now();	
  var elapsed = now - this.g_last;	
  this.g_last = now;
  this.g_stepCount = (this.g_stepCount +1)%1000;		// count 0,1,2,...999,0,1,2,...

  // Return the amount of time passed, in integer milliseconds
  return elapsed;
}

//Change Camera X Position 
var g_last2 = Date.now();
VBObox0.prototype.animate2 = function() {
//==============================================================================
  // Calculate the elapsed time
  var now = Date.now();
  var elapsed = now - g_last2;
  g_last2 = now;

  var newx = this.xCamPos + (this.moveSpeed*(Math.cos(Math.PI* (this.currentAngle/180)))) +
  (this.moveSpeed2*(Math.cos(Math.PI* ((this.currentAngle-90)/180))));// * elapsed);// / 1000.0;
  this.xCamPos = newx;
  return newx;
}
//Change Camera Y Position
var g_last3 = Date.now();
VBObox0.prototype.animate3 = function() {
//==============================================================================
  // Calculate the elapsed time
  var now = Date.now();
  var elapsed = now - g_last3;
  g_last3 = now;
console.log(this.moveSpeed);
 var newy = this.yCamPos + (this.moveSpeed*(Math.sin(Math.PI* (this.currentAngle/180)))) + 
     (this.moveSpeed2*(Math.sin(Math.PI* ((this.currentAngle-90)/180))));// * elapsed);// / 1000.0;
 	if(this.moveSpeed != 0) console.log(newy);

 //this.yCamPos = newy;
 //console.log('yo');
  return newy;
}

//Animate Angle of Camera
var g_lastA = Date.now();
VBObox0.prototype.animateA = function() {
//==============================================================================
  // Calculate the elapsed time
  var now = Date.now();
  var elapsed = now - g_lastA;
  g_lastA = now;


  var newAngle = this.currentAngle + ((this.ANGLE3_STEP * elapsed) / 1000.0);
  //console.log('ran =', newAngle);
  this.currentAngle = newAngle;
    //  	console.log(newAngle);
 //   console.log('hi');
  return this.currentAngle %= 360;

}*/
function makeGroundGrid() {
//==============================================================================
// Create a list of vertices that create a large grid of lines in the x,y plane
// centered at x=y=z=0.  Draw this shape using the GL_LINES primitive.

    var xymax = 1000.0;			// grid size; extends to cover +/-xymax in x and y.
    var xColr = new Float32Array([.4, .7, 0.1]);	// bright yellow
    var yColr = new Float32Array([0.1, 0, 0.8]);	// bright green.

    // Create an (global) array to hold this ground-plane's vertices:
    // draw a grid made of xcount+ycount lines; 2 vertices per line.

    var xgap = xymax / (xcount - 1);		// HALF-spacing between lines in x,y;
    var ygap = xymax / (ycount - 1);		// (why half? because v==(0line number/2))

    // First, step thru x values as we make vertical lines of constant-x:
    for (v = 0, j = 0; v < 2 * xcount; v++, j += floatsPerVertex) {
        if (v % 2 == 0) {	// put even-numbered vertices at (xnow, -xymax, 0)
            gndVerts[j] = -xymax + (v) * xgap;	// x
            gndVerts[j + 1] = -xymax;								// y
            gndVerts[j + 2] = -1.5;									// z
            gndVerts[j + 3] = 1.0;									// w.
        } else {				// put odd-numbered vertices at (xnow, +xymax, 0).
            gndVerts[j] = -xymax + (v - 1) * xgap;	// x
            gndVerts[j + 1] = xymax;								// y
            gndVerts[j + 2] = -1.5;									// z
            gndVerts[j + 3] = 1.0;									// w.
        }
        gndVerts[j + 4] = xColr[0];			// red
        gndVerts[j + 5] = xColr[1];			// grn
        gndVerts[j + 6] = xColr[2];			// blu
        /*		gndVerts[j+7] = 0;			// red
                gndVerts[j+8] = 0;			// grn
                gndVerts[j+9] = 1;*/
    }
    // Second, step thru y values as wqe make horizontal lines of constant-y:
    // (don't re-initialize j--we're adding more vertices to the array)
    for (v = 0; v < 2 * ycount; v++, j += floatsPerVertex) {
        if (v % 2 == 0) {		// put even-numbered vertices at (-xymax, ynow, 0)
            gndVerts[j] = -xymax;								// x
            gndVerts[j + 1] = -xymax + (v) * ygap;	// y
            gndVerts[j + 2] = -1.5;									// z
            gndVerts[j + 3] = 1.0;									// w.
        } else {					// put odd-numbered vertices at (+xymax, ynow, 0).
            gndVerts[j] = xymax;								// x
            gndVerts[j + 1] = -xymax + (v - 1) * ygap;	// y
            gndVerts[j + 2] = -1.5;									// z
            gndVerts[j + 3] = 1.0;									// w.
        }
        gndVerts[j + 4] = yColr[0];			// red
        gndVerts[j + 5] = yColr[1];			// grn
        gndVerts[j + 6] = yColr[2];			// blu
        /*		gndVerts[j+7] = 0;			// red
                gndVerts[j+8] = 0;			// grn
                gndVerts[j+9] = 1;*/
    }

    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;

    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;

    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = 0;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = axisSize;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = 0;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
    gndVerts[j] = axisSize;								// x
    gndVerts[j + 1] = 0;	// y
    gndVerts[j + 2] = axisSize;									// z
    gndVerts[j + 3] = 1.0;
    gndVerts[j + 4] = yColr[0];			// red
    gndVerts[j + 5] = xColr[1];			// grn
    gndVerts[j + 6] = yColr[2];			// blu
    j += floatsPerVertex;
}