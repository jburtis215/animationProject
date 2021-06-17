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
function VBObox3() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox' object  that holds all data and fcns 
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate set of shaders.
    this.VERT_SRC =	//--------------------- VERTEX SHADER source code
        'precision mediump float;\n' +				// req'd in OpenGL ES if we use 'float'
        //
        'uniform   int u_runMode1; \n' +
        'uniform float maxDensity; \n' +
        'uniform float minDensity; \n' +
        'uniform   mat4 u_MvpMat1; \n' +    //camera position
        'attribute vec3 a_Pos1;\n' +
   //     'attribute vec3 a_Colr1;\n' +
        'varying   vec4 v_Colr1; \n' +
        'varying   float r;\n' +
        'varying   float g;\n' +
        'varying   float b;\n' +
        'varying   float ratio;\n'  +
        'attribute float a_Density1; \n' +
        'attribute float a_diam1; \n' +
        'void main() {\n' +
        '  gl_PointSize = a_diam1;\n' +
        '	 gl_Position = u_MvpMat1 * (vec4(a_Pos1, 1)); \n' +// + u_ballShift1)
    //    'v_Colr1 = vec4(a_Colr1, 1.0); \n' +	// green: >3==run
        '   ratio = 2.0 * (a_Density1-minDensity) / (maxDensity - minDensity); \n' +
        '    b = max(0.0, (1.0 - ratio)); \n' +
        '    r = max(0.0, (ratio - 1.0)); \n' +
        '    g = 1.0- b - r; \n' +
        '   v_Colr1 = vec4(r, g, b, 1.0);\n' +

        '} \n';

    this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code
        'precision mediump float;\n' +
        'varying vec4 v_Colr1; \n' +
        'void main() {\n' +
        '  float dist = distance(gl_PointCoord, vec2(0.5, 0.5)); \n' +
        //  '  if(dist < 0.5) { \n' +
        '  	gl_FragColor = vec4((1.0-2.0*dist)*v_Colr1.rgb, 1.0);\n' +
        //	'  } else { discard; }\n' +
        '}\n';


    this.vboLoc;										// Vertex Buffer Object location# on the GPU
    this.vertices;
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


VBObox3.prototype.init = function (myGL) {
//=============================================================================
// Create, compile, link this VBObox object's shaders to an executable 'program'
// ready for use in the GPU.  Create and fill a Float32Array that holds all VBO 
// vertices' values; create a new VBO on the GPU and fill it with those values. 
// Find the GPU location of	all our shaders' attribute- and uniform-variables; 
// assign the correct portions of VBO contents as the data source for each 
// attribute, and transfer current values to the GPU for each uniform variable.
// (usually called only once, within main()) 

//partSys.s1 = new Float32Array(partSys.partCount*PART_MAXVAR);
    this.FSIZE = partSys1.s0.BYTES_PER_ELEMENT;	// memory needed to store an s0 array element.

// this.PartSys_init();
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
        partSys3.s0, 		// JavaScript Float32Array
        gl.DYNAMIC_DRAW);			// Usage hint.

// Find & Set All Attributes:------------------------------
    // a) Get the GPU location for each attribute var used in our shaders:
    this.a_PosLoc = gl.getAttribLocation(this.shaderLoc, 'a_Pos1');
    if (this.a_PosLoc < 0) {
        console.log(this.constructor.name +
            '.init() Failed to get GPU location of attribute a_Pos1');
        return -1;	// error exit.
    }

 /*   this.a_ColrLoc = myGL.getAttribLocation(this.shaderLoc, 'a_Colr1');
    if (this.a_ColrLoc < 0) {
        console.log(this.constructor.name +
            '.init() failed to get the GPU location of attribute a_Colr1');
        return -1;	// error exit.
    }*/
    this.a_diamID = myGL.getAttribLocation(this.shaderLoc, 'a_diam1');
    if (this.a_diamID < 0) {
        console.log(this.constructor.name +
            '.init() failed to get the GPU location of attribute a_diam1');
        return -1;	// error exit.
    }
    this.a_Density = myGL.getAttribLocation(this.shaderLoc, 'a_Density1');
    if (this.a_Density < 0) {
        console.log(this.constructor.name +
            '.init() failed to get the GPU location of attribute a_Density1');
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
  //  gl.vertexAttribPointer(this.a_ColrLoc, 3, gl.FLOAT, false,
    //    PART_MAXVAR * this.FSIZE, PART_R * this.FSIZE);
    myGL.vertexAttribPointer(this.a_diamID, 1, gl.FLOAT, false,
        PART_MAXVAR * this.FSIZE, PART_DIAM * this.FSIZE);
    // c) Enable this assignment of the attribute to its' VBO source:
    myGL.enableVertexAttribArray(this.a_PosLoc);
 //   myGL.enableVertexAttribArray(this.a_ColrLoc);
    myGL.enableVertexAttribArray(this.a_diamID);
// Find All Uniforms:--------------------------------
//Get GPU storage location for each uniform var used in our shader programs: 
    this.u_MvpMatLoc = myGL.getUniformLocation(this.shaderLoc, 'u_MvpMat1');
    if (!this.u_MvpMatLoc) {
        console.log(this.constructor.name +
            '.init() failed to get GPU location for u_MvpMat1 uniform');
        return;
    }

    this.maxDensity = myGL.getUniformLocation(this.shaderLoc, 'maxDensity');
    if (!this.maxDensity) {
        console.log(this.constructor.name +
            '.init() failed to get GPU location for maxDensity uniform');
        return;
    }
    this.minDensity = myGL.getUniformLocation(this.shaderLoc, 'minDensity');
    if (!this.minDensity) {
        console.log(this.constructor.name +
            '.init() failed to get GPU location for minDensity uniform');
        return;
    }
}


VBObox3.prototype.adjust1 = function (mygl, partSys) {
    mygl.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    partSys.minDensity = 1000000;
    partSys.maxDensity = 1;
    for (var i = 0; i < partSys.partCount; i++) {
        this.applyForces(mygl, i, partSys);
        this.dotFinder(mygl, i, partSys);
        this.render(mygl, i, partSys);
        this.PartSys_constrain1(i, partSys);
    }
    partSys.buffer = partSys.s0.slice();
    partSys.s0 = partSys.s1.slice();
    partSys.s1 = partSys.buffer.slice();
    this.vertices = partSys.s0;
    this.PartSys_render(gl, partSys);

}

VBObox3.prototype.applyForces = function (mygl, i, partSys) {

    var rad = 0;
    var pOff = i * PART_MAXVAR;
    partSys.s0[pOff + PART_X_FTOT] = 0;
    partSys.s0[pOff + PART_Y_FTOT] = 0;
    partSys.s0[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_E](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_P](pOff, partSys.s0);
    partSys.cforcer[F_WIND](pOff, partSys.s0);
    partSys.cforcer[F_BUBBLE](pOff, partSys.s0);
    partSys.cforcer[F_DRAG](pOff, partSys.s0);
    partSys.cforcer[F_SPRING](pOff, partSys.s0);
    partSys.cforcer[F_SPH](pOff, partSys.s0);
    partSys.cforcer[F_CHARGE](pOff, partSys.s0);
    partSys.cforcer[F_FIRE](pOff, partSys.s0);
   // console.log("this is the new density" + partSys.maxDensity);
}

VBObox3.prototype.dotFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.sdot[pOff + PART_XPOS] = partSys.s0[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
    partSys.sdot[pOff + PART_YPOS] = partSys.s0[pOff + PART_YVEL];
    partSys.sdot[pOff + PART_ZPOS] = partSys.s0[pOff + PART_ZVEL];
    partSys.sdot[pOff + PART_XVEL] = partSys.s0[pOff + PART_X_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.sdot[pOff + PART_YVEL] = partSys.s0[pOff + PART_Y_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.sdot[pOff + PART_ZVEL] = partSys.s0[pOff + PART_Z_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.sdot[pOff + PART_X_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Y_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Z_FTOT] = 0.0;
    partSys.sdot[pOff + PART_R] = 0;
    partSys.sdot[pOff + PART_G] = partSys.s0[pOff + PART_G_VEL];
    partSys.sdot[pOff + PART_B] = partSys.s0[pOff + PART_B_VEL];
    partSys.sdot[pOff + PART_MASS] = partSys.sdot[pOff + PART_MASS_VEL];
    partSys.sdot[pOff + PART_DIAM] = 0;
    partSys.sdot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.sdot[pOff + PART_AGE] = .08;
    //		partSys.sdot[pOff + PART_AGE]      = (i/partSys.partCount) * 10;  // # of frame-times since creation/initialization
    partSys.sdot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.sdot[pOff + PART_MASS_VEL] = partSys.s0[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
    partSys.sdot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.sdot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.sdot[pOff + PART_G_VEL] = partSys.s0[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
    partSys.sdot[pOff + PART_B_VEL] = partSys.s0[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
    partSys.sdot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.sdot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.sdot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
    partSys.sdot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
}

VBObox3.prototype.render = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s1[pOff + PART_XPOS] = (partSys.sdot[pOff + PART_XPOS] / timeStep) + partSys.s0[pOff + PART_XPOS];
    partSys.s1[pOff + PART_YPOS] = (partSys.sdot[pOff + PART_YPOS] / timeStep) + partSys.s0[pOff + PART_YPOS];
    partSys.s1[pOff + PART_ZPOS] = (partSys.sdot[pOff + PART_ZPOS] / timeStep) + partSys.s0[pOff + PART_ZPOS];
    partSys.s1[pOff + PART_XVEL] = (partSys.sdot[pOff + PART_XVEL] / timeStep) + partSys.s0[pOff + PART_XVEL];
    partSys.s1[pOff + PART_YVEL] = (partSys.sdot[pOff + PART_YVEL] / timeStep) + partSys.s0[pOff + PART_YVEL];
    partSys.s1[pOff + PART_ZVEL] = (partSys.sdot[pOff + PART_ZVEL] / timeStep) + partSys.s0[pOff + PART_ZVEL];
    partSys.s1[pOff + PART_R] = (partSys.sdot[pOff + PART_R] / timeStep) + partSys.s0[pOff + PART_R];
    partSys.s1[pOff + PART_G] = (partSys.sdot[pOff + PART_G] / timeStep) + partSys.s0[pOff + PART_G];
    partSys.s1[pOff + PART_B] = (partSys.sdot[pOff + PART_B] / timeStep) + partSys.s0[pOff + PART_B];
    partSys.s1[pOff + PART_R_VEL] = (partSys.sdot[pOff + PART_R_VEL] / timeStep) + partSys.s0[pOff + PART_R_VEL];
    partSys.s1[pOff + PART_G_VEL] = (partSys.sdot[pOff + PART_G_VEL] / timeStep) + partSys.s0[pOff + PART_G_VEL];
    partSys.s1[pOff + PART_B_VEL] = (partSys.sdot[pOff + PART_B_VEL] / timeStep) + partSys.s0[pOff + PART_B_VEL];
    partSys.s1[pOff + PART_MASS] = (partSys.sdot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS];
    partSys.s1[pOff + PART_AGE] = partSys.sdot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];
    partSys.s1[pOff + PART_DENSITY] = findDensity(partSys.s1, pOff, partSys.partCount);
  //  partSys.s1[pOff + PART_PRESSURE] =
    if (partSys.maxDensity < partSys.s1[pOff + PART_DENSITY]){
        partSys.maxDensity = partSys.s1[pOff + PART_DENSITY];
    }
    if (partSys.minDensity > partSys.s1[pOff + PART_DENSITY]){
        partSys.minDensity = partSys.s1[pOff + PART_DENSITY];
    }
}

VBObox3.prototype.PartSys_constrain1 = function (i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.c0[0](i, partSys.s1, partSys.s0);
    partSys.c0[1](pOff, partSys.s1);
}

VBObox3.prototype.PartSys_render = function (myGL, partSys) {

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
    // 	myGL.uniform1i(this.u_runModeID, myRunMode);		// run/step/pause the particle system

    myGL.uniform1f(this.maxDensity,
        partSys.maxDensity);
    myGL.uniform1f(this.minDensity,
        partSys.minDensity);
    myGL.bindBuffer(myGL.ARRAY_BUFFER,	// GLenum 'target' for this GPU buffer
        this.vboLoc);			// the ID# the GPU uses for this buffer.
    // (Here's how to use the almost-identical OpenGL version of this function:
    //		http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribPointer.xml )
    // b) Re-connect data paths from VBO to each shader attribute:
    myGL.vertexAttribPointer(this.a_PosLoc, 3, myGL.FLOAT, false,
        PART_MAXVAR * this.FSIZE, 0);		// stride, offset
   // myGL.vertexAttribPointer(this.a_ColrLoc, 3, myGL.FLOAT, false,
     //   PART_MAXVAR * this.FSIZE, PART_R * this.FSIZE); // stride, offset
    myGL.vertexAttribPointer(this.a_diamID, 1, myGL.FLOAT, false,
        PART_MAXVAR * this.FSIZE, PART_DIAM * this.FSIZE); // stride, offset

    myGL.vertexAttribPointer(this.a_Density, 1, myGL.FLOAT, false,
        PART_MAXVAR * this.FSIZE, PART_DENSITY * this.FSIZE); // stride, offset
    // c) enable the newly-re-assigned attributes:
    myGL.enableVertexAttribArray(this.a_PosLoc);
 //   myGL.enableVertexAttribArray(this.a_ColrLoc);
    myGL.enableVertexAttribArray(this.a_diamID);
    myGL.enableVertexAttribArray(this.a_Density);

    // ----------------------------Draw the contents of the currently-bound VBO:
    myGL.bufferSubData(myGL.ARRAY_BUFFER, 0, partSys.s0);
    myGL.drawArrays(myGL.POINTS, 0, partSys.partCount);

}


VBObox3.prototype.reset = function (partSys) {
    /*for (i = 0; i < partSys.partCount; i++){
                if(s0[i* PART_MAXVAR + PART_XVEL] > 0.0) s0[i* PART_MAXVAR + PART_XVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_XVEL] -= INIT_VEL;
                if(s0[i* PART_MAXVAR + PART_YVEL] > 0.0) s0[i* PART_MAXVAR + PART_YVEL] += INIT_VEL*.8; else s0[i* PART_MAXVAR + PART_YVEL] -= INIT_VEL*.8;
                if(s0[i* PART_MAXVAR + PART_ZVEL] > 0.0) s0[i* PART_MAXVAR + PART_ZVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_ZVEL] -= INIT_VEL;
                }*/
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s0[pOff + PART_XVEL] > 0) {
            partSys.s0[pOff + PART_XVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_XVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_YVEL] > 0) {
            partSys.s0[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_ZVEL] > 0) {
            partSys.s0[pOff + PART_ZVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_ZVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
    }
}
VBObox3.prototype.resetWind = function (partSys) {
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s0[pOff + PART_XVEL] > 0) {
            partSys.s0[pOff + PART_XVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_XVEL] -= (0.3 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_YVEL] > 0) {
            partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_ZVEL] > 0) {
            partSys.s0[pOff + PART_ZVEL] += (0.8 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_ZVEL] -= (0.8 + 0.8 * Math.random()) * INIT_VEL;
    }
}
VBObox3.prototype.resetSpring = function (partSys) {
    for (var i = 1; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        partSys.s0[pOff + PART_ZPOS] -= springRun * (.002 * Math.pow(i, 2));
    }
}

function VBObox4() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox' object  that holds all data and fcns 
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate set of shaders.
    this.VERT_SRC =	//--------------------- VERTEX SHADER source code
        'precision mediump float;\n' +				// req'd in OpenGL ES if we use 'float'
        //
        'uniform   int u_runMode1; \n' +					// particle system state:
        // 0=reset; 1= pause; 2=step; 3=run
        // 'uniform	 vec4 u_ballShift1; \n' +			// single bouncy-ball's movement
        'uniform   mat4 u_MvpMat1; \n' +    //camera position
        'attribute vec3 a_Pos1;\n' +
        'attribute vec3 a_Colr1;\n' +
        'varying   vec4 v_Colr1; \n' +
        'attribute float a_diam1; \n' +
        'void main() {\n' +
        '  gl_PointSize = a_diam1;\n' +
        '	 gl_Position = u_MvpMat1 * (vec4(a_Pos1, 1)); \n' +// + u_ballShift1)

        'v_Colr1 = vec4(a_Colr1, 1.0); \n' +	// green: >3==run
        '} \n';

    this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code
        'precision mediump float;\n' +
        'varying vec4 v_Colr1; \n' +
        'void main() {\n' +
        '  float dist = distance(gl_PointCoord, vec2(0.5, 0.5)); \n' +
        //  '  if(dist < 0.5) { \n' +
        '  	gl_FragColor = vec4((1.0-2.0*dist)*v_Colr1.rgb, 1.0);\n' +
        //	'  } else { discard; }\n' +
        '}\n';


    this.vboLoc;										// Vertex Buffer Object location# on the GPU
    this.vertices;
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


VBObox4.prototype.init = function (myGL) {
//=============================================================================
// Create, compile, link this VBObox object's shaders to an executable 'program'
// ready for use in the GPU.  Create and fill a Float32Array that holds all VBO 
// vertices' values; create a new VBO on the GPU and fill it with those values. 
// Find the GPU location of	all our shaders' attribute- and uniform-variables; 
// assign the correct portions of VBO contents as the data source for each 
// attribute, and transfer current values to the GPU for each uniform variable.
// (usually called only once, within main()) 

//partSys.s1 = new Float32Array(partSys.partCount*PART_MAXVAR);
    this.FSIZE = partSys1.s0.BYTES_PER_ELEMENT;	// memory needed to store an s0 array element.

// this.PartSys_init();
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
        partSys1.s0, 		// JavaScript Float32Array
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

}


VBObox4.prototype.adjust1 = function (mygl, partSys) {
    mygl.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    for (var i = 0; i < partSys.partCount; i++) {
        this.applyForces(mygl, i, partSys);
        this.halfdotFinder(mygl, i, partSys);
        this.midPointFinder(mygl, i, partSys);
        this.applyMidForces(mygl, i, partSys);
        this.midDotFinder(mygl, i, partSys);
        this.render(mygl, i, partSys);
        this.PartSys_constrain1(i, partSys);
    }
    partSys.buffer = partSys.s0.slice();
    partSys.s0 = partSys.s1.slice();
    partSys.s1 = partSys.buffer.slice();
    this.vertices = partSys.s0;
    this.PartSys_render(gl, partSys);

}

VBObox4.prototype.applyForces = function (mygl, i, partSys) {

    var rad = 0;
    var pOff = i * PART_MAXVAR;
    partSys.s0[pOff + PART_X_FTOT] = 0;
    partSys.s0[pOff + PART_Y_FTOT] = 0;
    partSys.s0[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_E](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_P](pOff, partSys.s0);
    partSys.cforcer[F_WIND](pOff, partSys.s0);
    partSys.cforcer[F_BUBBLE](pOff, partSys.s0);
    partSys.cforcer[F_DRAG](pOff, partSys.s0);
    partSys.cforcer[F_SPRING](pOff, partSys.s0);
    partSys.cforcer[F_SPH](pOff, partSys.s0);
    partSys.cforcer[F_CHARGE](pOff, partSys.s0);
    partSys.cforcer[F_FIRE](pOff, partSys.s0);

}

VBObox4.prototype.halfdotFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.sdot[pOff + PART_XPOS] = partSys.s0[pOff + PART_XVEL] / 2;		// 0.0 <= randomRound() < 1.0
    partSys.sdot[pOff + PART_YPOS] = partSys.s0[pOff + PART_YVEL] / 2;
    partSys.sdot[pOff + PART_ZPOS] = partSys.s0[pOff + PART_ZVEL] / 2;
    partSys.sdot[pOff + PART_XVEL] = partSys.s0[pOff + PART_X_FTOT] * partSys.s1[pOff + PART_MASS] / 2;
    partSys.sdot[pOff + PART_YVEL] = partSys.s0[pOff + PART_Y_FTOT] * partSys.s1[pOff + PART_MASS] / 2;
    partSys.sdot[pOff + PART_ZVEL] = partSys.s0[pOff + PART_Z_FTOT] * partSys.s1[pOff + PART_MASS] / 2;
    partSys.sdot[pOff + PART_X_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Y_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Z_FTOT] = 0.0;
    partSys.sdot[pOff + PART_R] = 0;
    partSys.sdot[pOff + PART_G] = partSys.s0[pOff + PART_G_VEL] / 2;
    partSys.sdot[pOff + PART_B] = partSys.s0[pOff + PART_B_VEL] / 2;
    partSys.sdot[pOff + PART_MASS] = partSys.s0[pOff + PART_MASS_VEL] / 2;
    partSys.sdot[pOff + PART_DIAM] = 0;
    partSys.sdot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.sdot[pOff + PART_AGE] = .08 / 2;
    partSys.sdot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.sdot[pOff + PART_MASS_VEL] = partSys.s0[pOff + PART_MASS_FTOT] / 2;  // time-rate-of-change of mass.
    partSys.sdot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.sdot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.sdot[pOff + PART_G_VEL] = partSys.s0[pOff + PART_G_FTOT] / 2;  // time-rate-of-change of color:grn
    partSys.sdot[pOff + PART_B_VEL] = partSys.s0[pOff + PART_B_FTOT] / 2; // time-rate-of-change of color:blu
    partSys.sdot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.sdot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.sdot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu

}

VBObox4.prototype.midPointFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.sM[pOff + PART_XPOS] = (partSys.sdot[pOff + PART_XPOS] / timeStep) + partSys.s0[pOff + PART_XPOS];
    partSys.sM[pOff + PART_YPOS] = (partSys.sdot[pOff + PART_YPOS] / timeStep) + partSys.s0[pOff + PART_YPOS];
    partSys.sM[pOff + PART_ZPOS] = (partSys.sdot[pOff + PART_ZPOS] / timeStep) + partSys.s0[pOff + PART_ZPOS];
    partSys.sM[pOff + PART_XVEL] = (partSys.sdot[pOff + PART_XVEL] / timeStep) + partSys.s0[pOff + PART_XVEL];
    partSys.sM[pOff + PART_YVEL] = (partSys.sdot[pOff + PART_YVEL] / timeStep) + partSys.s0[pOff + PART_YVEL];
    partSys.sM[pOff + PART_ZVEL] = (partSys.sdot[pOff + PART_ZVEL] / timeStep) + partSys.s0[pOff + PART_ZVEL];
    partSys.sM[pOff + PART_R] = (partSys.sdot[pOff + PART_R] / timeStep) + partSys.s0[pOff + PART_R];
    partSys.sM[pOff + PART_G] = (partSys.sdot[pOff + PART_G] / timeStep) + partSys.s0[pOff + PART_G];
    partSys.sM[pOff + PART_B] = (partSys.sdot[pOff + PART_B] / timeStep) + partSys.s0[pOff + PART_B];
    partSys.sM[pOff + PART_R_VEL] = (partSys.sdot[pOff + PART_R_VEL] / timeStep) + partSys.s0[pOff + PART_R_VEL];
    partSys.sM[pOff + PART_G_VEL] = (partSys.sdot[pOff + PART_G_VEL] / timeStep) + partSys.s0[pOff + PART_G_VEL];
    partSys.sM[pOff + PART_B_VEL] = (partSys.sdot[pOff + PART_B_VEL] / timeStep) + partSys.s0[pOff + PART_B_VEL];
    partSys.sM[pOff + PART_MASS] = (partSys.sdot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS];
    partSys.sM[pOff + PART_AGE] = partSys.sdot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];

}

VBObox4.prototype.applyMidForces = function (mygl, i, partSys) {

    var rad = 0;

    var pOff = i * PART_MAXVAR;
    partSys.sM[pOff + PART_X_FTOT] = 0;
    partSys.sM[pOff + PART_Y_FTOT] = 0;
    partSys.sM[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.sM);
    partSys.cforcer[F_GRAV_E](pOff, partSys.sM);
    partSys.cforcer[F_GRAV_P](pOff, partSys.sM);
    partSys.cforcer[F_WIND](pOff, partSys.sM);
    partSys.cforcer[F_BUBBLE](pOff, partSys.sM);
    partSys.cforcer[F_DRAG](pOff, partSys.sM);

    partSys.cforcer[F_SPRING](pOff, partSys.sM);
    partSys.cforcer[F_SPH](pOff, partSys.sM);
    partSys.cforcer[F_CHARGE](pOff, partSys.sM);
    partSys.cforcer[F_FIRE](pOff, partSys.sM);

}

VBObox4.prototype.midDotFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s2dot[pOff + PART_XPOS] = partSys.sM[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
    partSys.s2dot[pOff + PART_YPOS] = partSys.sM[pOff + PART_YVEL];
    partSys.s2dot[pOff + PART_ZPOS] = partSys.sM[pOff + PART_ZVEL];
    partSys.s2dot[pOff + PART_XVEL] = partSys.sM[pOff + PART_X_FTOT] * partSys.sM[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_YVEL] = partSys.sM[pOff + PART_Y_FTOT] * partSys.sM[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_ZVEL] = partSys.sM[pOff + PART_Z_FTOT] * partSys.sM[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_X_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_Y_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_Z_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_R] = 0;
    partSys.s2dot[pOff + PART_G] = partSys.sM[pOff + PART_G_VEL];
    partSys.s2dot[pOff + PART_B] = partSys.sM[pOff + PART_B_VEL];
    partSys.s2dot[pOff + PART_MASS] = partSys.sM[pOff + PART_MASS_VEL];
    partSys.s2dot[pOff + PART_DIAM] = 0;
    partSys.s2dot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.s2dot[pOff + PART_AGE] = .08;
    //		partSys.s2dot[pOff + PART_AGE]      = (i/partSys.partCount) * 10;  // # of frame-times since creation/initialization
    partSys.s2dot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.s2dot[pOff + PART_MASS_VEL] = partSys.sM[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
    partSys.s2dot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.s2dot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.s2dot[pOff + PART_G_VEL] = partSys.sM[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
    partSys.s2dot[pOff + PART_B_VEL] = partSys.sM[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
    partSys.s2dot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.s2dot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.s2dot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu

}
VBObox4.prototype.render = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s1[pOff + PART_XPOS] = (partSys.s2dot[pOff + PART_XPOS] / timeStep) + partSys.s0[pOff + PART_XPOS];
    partSys.s1[pOff + PART_YPOS] = (partSys.s2dot[pOff + PART_YPOS] / timeStep) + partSys.s0[pOff + PART_YPOS];
    partSys.s1[pOff + PART_ZPOS] = (partSys.s2dot[pOff + PART_ZPOS] / timeStep) + partSys.s0[pOff + PART_ZPOS];
    partSys.s1[pOff + PART_XVEL] = (partSys.s2dot[pOff + PART_XVEL] / timeStep) + partSys.s0[pOff + PART_XVEL];
    partSys.s1[pOff + PART_YVEL] = (partSys.s2dot[pOff + PART_YVEL] / timeStep) + partSys.s0[pOff + PART_YVEL];
    partSys.s1[pOff + PART_ZVEL] = (partSys.s2dot[pOff + PART_ZVEL] / timeStep) + partSys.s0[pOff + PART_ZVEL];
    partSys.s1[pOff + PART_R] = (partSys.s2dot[pOff + PART_R] / timeStep) + partSys.s0[pOff + PART_R];
    partSys.s1[pOff + PART_G] = (partSys.s2dot[pOff + PART_G] / timeStep) + partSys.s0[pOff + PART_G];
    partSys.s1[pOff + PART_B] = (partSys.s2dot[pOff + PART_B] / timeStep) + partSys.s0[pOff + PART_B];
    partSys.s1[pOff + PART_R_VEL] = (partSys.s2dot[pOff + PART_R_VEL] / timeStep) + partSys.s0[pOff + PART_R_VEL];
    partSys.s1[pOff + PART_G_VEL] = (partSys.s2dot[pOff + PART_G_VEL] / timeStep) + partSys.s0[pOff + PART_G_VEL];
    partSys.s1[pOff + PART_B_VEL] = (partSys.s2dot[pOff + PART_B_VEL] / timeStep) + partSys.s0[pOff + PART_B_VEL];
    partSys.s1[pOff + PART_MASS] = (partSys.s2dot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS];
    partSys.s1[pOff + PART_AGE] = partSys.s2dot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];

}

VBObox4.prototype.PartSys_constrain1 = function (i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.c0[0](i, partSys.s1, partSys.s0);
    partSys.c0[1](pOff, partSys.s1);
}

VBObox4.prototype.PartSys_render = function (myGL, partSys) {

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

    // 	myGL.uniform1i(this.u_runModeID, myRunMode);		// run/step/pause the particle system

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
    myGL.bufferSubData(myGL.ARRAY_BUFFER, 0, partSys.s0);
    myGL.drawArrays(myGL.POINTS, 0, partSys.partCount);

}


VBObox4.prototype.reset = function (partSys) {
    /*for (i = 0; i < partSys.partCount; i++){
                if(s0[i* PART_MAXVAR + PART_XVEL] > 0.0) s0[i* PART_MAXVAR + PART_XVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_XVEL] -= INIT_VEL;
                if(s0[i* PART_MAXVAR + PART_YVEL] > 0.0) s0[i* PART_MAXVAR + PART_YVEL] += INIT_VEL*.8; else s0[i* PART_MAXVAR + PART_YVEL] -= INIT_VEL*.8;
                if(s0[i* PART_MAXVAR + PART_ZVEL] > 0.0) s0[i* PART_MAXVAR + PART_ZVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_ZVEL] -= INIT_VEL;
                }*/
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s0[pOff + PART_XVEL] > 0) {
            partSys.s0[pOff + PART_XVEL] += (-0.4 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_XVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_YVEL] > 0) {
            partSys.s0[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_ZVEL] > 0) {
            partSys.s0[pOff + PART_ZVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_ZVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
    }
}

VBObox4.prototype.resetWind = function (partSys) {
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s0[pOff + PART_XVEL] > 0) {
            partSys.s0[pOff + PART_XVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_XVEL] -= (0.3 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_YVEL] > 0) {
            partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
    }
}

VBObox4.prototype.resetSpring = function (partSys) {
    for (var i = 1; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        partSys.s0[pOff + PART_ZPOS] -= springRun * (.002 * Math.pow(i, 2));
        /*			}
                else partSys.s0[pOff + PART_XVEL] -= (0.3 + 0.8*Math.random())*INIT_VEL;
                if(  partSys.s0[pOff + PART_YVEL] > 0) {
                       partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8*Math.random())*INIT_VEL;
                    }
                else partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8*Math.random())*INIT_VEL;
            }*/
    }
}

function VBObox5() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox' object  that holds all data and fcns 
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate set of shaders.
    this.VERT_SRC =	//--------------------- VERTEX SHADER source code
        'precision mediump float;\n' +				// req'd in OpenGL ES if we use 'float'
        //
        'uniform   int u_runMode1; \n' +					// particle system state:
        // 0=reset; 1= pause; 2=step; 3=run
        // 'uniform	 vec4 u_ballShift1; \n' +			// single bouncy-ball's movement
        'uniform   mat4 u_MvpMat1; \n' +    //camera position
        'attribute vec3 a_Pos1;\n' +
        'attribute vec3 a_Colr1;\n' +
        'varying   vec4 v_Colr1; \n' +
        'attribute float a_diam1; \n' +
        'void main() {\n' +
        '  gl_PointSize = a_diam1;\n' +
        '	 gl_Position = u_MvpMat1 * (vec4(a_Pos1, 1)); \n' +// + u_ballShift1)

        // Let u_runMode determine particle color:
        //  '  if(u_runMode1 == 0) { \n' +
        //	'	   v_Colr1 = vec4(1.0, 0.0, 0.0, 1.0);	\n' +		// red: 0==reset
        //	'  	 } \n' +
        //	'  else if(u_runMode1 == 1) {  \n' +
        //	'    v_Colr1 = vec4(1.0, 1.0, 0.0, 1.0); \n' +	// yellow: 1==pause
        //	'    }  \n' +
        //	'  else if(u_runMode1 == 2) { \n' +
        //	'    v_Colr1 = vec4(1.0, 1.0, 1.0, 1.0); \n' +	// white: 2==step
        //'    } \n' +
        //'  else { \n' +
        'v_Colr1 = vec4(a_Colr1, 1.0); \n' +	// green: >3==run
        //	'		 } \n' +
        '} \n';

    this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code
        'precision mediump float;\n' +
        'varying vec4 v_Colr1; \n' +
        'void main() {\n' +
        '  float dist = distance(gl_PointCoord, vec2(0.5, 0.5)); \n' +
        //  '  if(dist < 0.5) { \n' +
        '  	gl_FragColor = vec4((1.0-2.0*dist)*v_Colr1.rgb, 1.0);\n' +
        //	'  } else { discard; }\n' +
        '}\n';


    this.vboLoc;										// Vertex Buffer Object location# on the GPU
    this.vertices;
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


VBObox5.prototype.init = function (myGL) {
//=============================================================================
// Create, compile, link this VBObox object's shaders to an executable 'program'
// ready for use in the GPU.  Create and fill a Float32Array that holds all VBO 
// vertices' values; create a new VBO on the GPU and fill it with those values. 
// Find the GPU location of	all our shaders' attribute- and uniform-variables; 
// assign the correct portions of VBO contents as the data source for each 
// attribute, and transfer current values to the GPU for each uniform variable.
// (usually called only once, within main()) 

//partSys.s1 = new Float32Array(partSys.partCount*PART_MAXVAR);
    this.FSIZE = partSys1.s0.BYTES_PER_ELEMENT;	// memory needed to store an s0 array element.

// this.PartSys_init();
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
        partSys1.s0, 		// JavaScript Float32Array
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


}


VBObox5.prototype.adjust1 = function (mygl, partSys) {
    mygl.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    for (var i = 0; i < partSys.partCount; i++) {
        this.applyForces(mygl, i, partSys);
        this.dotFinder(mygl, i, partSys);
        this.computer(mygl, i, partSys);
        this.apply2Forces(mygl, i, partSys);
        this.Dot2Finder(mygl, i, partSys);
        this.veloComputer(mygl, i, partSys);
        this.PartSys_constrain1(i, partSys);
    }
    partSys.buffer = partSys.s0.slice();
    partSys.s0 = partSys.s1.slice();
    partSys.s1 = partSys.buffer.slice();
    this.vertices = partSys.s0;
    this.PartSys_render(gl, partSys);

}

/*VBObox5.prototype.oldPointFinder = function(mygl, i, partSys){
				var pOff = i*PART_MAXVAR;
				partSys.s0[pOff + PART_XPOS] = partSys.s1[pOff + PART_XPOS] - (partSys.sdot[pOff + PART_XPOS] / timeStep);
				partSys.s0[pOff + PART_YPOS] = partSys.s1[pOff + PART_YPOS] - (partSys.sdot[pOff + PART_YPOS] / timeStep);
				partSys.s0[pOff + PART_ZPOS] = partSys.s1[pOff + PART_ZPOS] - (partSys.sdot[pOff + PART_ZPOS] / timeStep);
				partSys.s0[pOff + PART_XVEL] = partSys.s1[pOff + PART_XVEL] - (partSys.sdot[pOff + PART_XVEL] / timeStep);
				partSys.s0[pOff + PART_YVEL] = partSys.s1[pOff + PART_YVEL] - (partSys.sdot[pOff + PART_YVEL] / timeStep);
				partSys.s0[pOff + PART_ZVEL] = partSys.s1[pOff + PART_ZVEL] - (partSys.sdot[pOff + PART_ZVEL] / timeStep);
				partSys.s0[pOff + PART_R] =  partSys.s1[pOff + PART_R] - (partSys.sdot[pOff + PART_R] / timeStep);
				partSys.s0[pOff + PART_G] = partSys.s1[pOff + PART_G] - (partSys.sdot[pOff + PART_G] / timeStep);
				partSys.s0[pOff + PART_B] = partSys.s1[pOff + PART_B] - (partSys.sdot[pOff + PART_B] / timeStep);
				partSys.s0[pOff + PART_R_VEL] = partSys.s1[pOff + PART_R_VEL] - (partSys.sdot[pOff + PART_R_VEL] / timeStep);
				partSys.s0[pOff + PART_G_VEL] = partSys.s1[pOff + PART_G_VEL] - (partSys.sdot[pOff + PART_G_VEL] / timeStep);
				partSys.s0[pOff + PART_B_VEL] = partSys.s1[pOff + PART_B_VEL] - (partSys.sdot[pOff + PART_B_VEL] / timeStep);
				partSys.s0[pOff + PART_MASS] = partSys.s1[pOff + PART_MASS] - (partSys.sdot[pOff + PART_MASS] / timeStep);
				partSys.s0[pOff + PART_AGE] = partSys.s1[pOff + PART_AGE] - partSys.sdot[pOff + PART_AGE];
				
}*/
VBObox5.prototype.applyForces = function (mygl, i, partSys) {

    var rad = 0;
    var pOff = i * PART_MAXVAR;
    partSys.s0[pOff + PART_X_FTOT] = 0;
    partSys.s0[pOff + PART_Y_FTOT] = 0;
    partSys.s0[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_E](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_P](pOff, partSys.s0);
    partSys.cforcer[F_WIND](pOff, partSys.s0);
    partSys.cforcer[F_BUBBLE](pOff, partSys.s0);
    partSys.cforcer[F_DRAG](pOff, partSys.s0);
    partSys.cforcer[F_SPRING](pOff, partSys.s0);
    partSys.cforcer[F_SPH](pOff, partSys.s0);
    partSys.cforcer[F_CHARGE](pOff, partSys.s0);
    partSys.cforcer[F_FIRE](pOff, partSys.s0);

}

VBObox5.prototype.dotFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.sdot[pOff + PART_XPOS] = partSys.s0[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
    partSys.sdot[pOff + PART_YPOS] = partSys.s0[pOff + PART_YVEL];
    partSys.sdot[pOff + PART_ZPOS] = partSys.s0[pOff + PART_ZVEL];
    partSys.sdot[pOff + PART_XVEL] = partSys.s0[pOff + PART_X_FTOT] * partSys.s0[pOff + PART_MASS];
    partSys.sdot[pOff + PART_YVEL] = partSys.s0[pOff + PART_Y_FTOT] * partSys.s0[pOff + PART_MASS];
    partSys.sdot[pOff + PART_ZVEL] = partSys.s0[pOff + PART_Z_FTOT] * partSys.s0[pOff + PART_MASS];
    partSys.sdot[pOff + PART_X_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Y_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Z_FTOT] = 0.0;
    partSys.sdot[pOff + PART_R] = 0;
    partSys.sdot[pOff + PART_G] = partSys.s0[pOff + PART_G_VEL];
    partSys.sdot[pOff + PART_B] = partSys.s0[pOff + PART_B_VEL];
    partSys.sdot[pOff + PART_MASS] = partSys.sdot[pOff + PART_MASS_VEL];
    partSys.sdot[pOff + PART_DIAM] = 0;
    partSys.sdot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.sdot[pOff + PART_AGE] = .08;
    //		partSys.sdot[pOff + PART_AGE]      = (i/partSys.partCount) * 10;  // # of frame-times since creation/initialization
    partSys.sdot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.sdot[pOff + PART_MASS_VEL] = partSys.s0[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
    partSys.sdot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.sdot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.sdot[pOff + PART_G_VEL] = partSys.s0[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
    partSys.sdot[pOff + PART_B_VEL] = partSys.s0[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
    partSys.sdot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.sdot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.sdot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
}

VBObox5.prototype.computer = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    if (pOff == 0) console.log(partSys.s0[pOff + PART_MASS]);
    partSys.s1[pOff + PART_XPOS] = (partSys.sdot[pOff + PART_XPOS] / timeStep) + partSys.s0[pOff + PART_XPOS] + partSys.sdot[pOff + PART_XVEL] / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_YPOS] = (partSys.sdot[pOff + PART_YPOS] / timeStep) + partSys.s0[pOff + PART_YPOS] + partSys.sdot[pOff + PART_YVEL] / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_ZPOS] = (partSys.sdot[pOff + PART_ZPOS] / timeStep) + partSys.s0[pOff + PART_ZPOS] + partSys.sdot[pOff + PART_ZVEL] / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_R] = (partSys.sdot[pOff + PART_R] / timeStep) + partSys.s0[pOff + PART_R] + partSys.s0[pOff + PART_R_VEL] / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_G] = (partSys.sdot[pOff + PART_G] / timeStep) + partSys.s0[pOff + PART_G] + partSys.s0[pOff + PART_G_VEL] / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_B] = (partSys.sdot[pOff + PART_B] / timeStep) + partSys.s0[pOff + PART_B] + partSys.s0[pOff + PART_B_VEL] / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_MASS] = (partSys.s2dot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS];
    partSys.s1[pOff + PART_AGE] = partSys.s2dot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];
}

VBObox5.prototype.apply2Forces = function (mygl, i, partSys) {

    var rad = 0;
    var pOff = i * PART_MAXVAR;
    partSys.s1[pOff + PART_X_FTOT] = 0;
    partSys.s1[pOff + PART_Y_FTOT] = 0;
    partSys.s1[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.s1);
    partSys.cforcer[F_GRAV_E](pOff, partSys.s1);
    partSys.cforcer[F_GRAV_P](pOff, partSys.s1);
    partSys.cforcer[F_WIND](pOff, partSys.s1);
    partSys.cforcer[F_BUBBLE](pOff, partSys.s1);
    partSys.cforcer[F_DRAG](pOff, partSys.s1);
    partSys.cforcer[F_SPRING](pOff, partSys.s1);
    partSys.cforcer[F_SPH](pOff, partSys.s1);
    partSys.cforcer[F_CHARGE](pOff, partSys.s1);
    partSys.cforcer[F_FIRE](pOff, partSys.s1);

}

VBObox5.prototype.Dot2Finder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s2dot[pOff + PART_XPOS] = partSys.s1[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
    partSys.s2dot[pOff + PART_YPOS] = partSys.s1[pOff + PART_YVEL];
    partSys.s2dot[pOff + PART_ZPOS] = partSys.s1[pOff + PART_ZVEL];

    partSys.s2dot[pOff + PART_XVEL] = partSys.s1[pOff + PART_X_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_YVEL] = partSys.s1[pOff + PART_Y_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_ZVEL] = partSys.s1[pOff + PART_Z_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_X_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_Y_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_Z_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_R] = 0;
    partSys.s2dot[pOff + PART_G] = partSys.s1[pOff + PART_G_VEL];
    partSys.s2dot[pOff + PART_B] = partSys.s1[pOff + PART_B_VEL];
    partSys.s2dot[pOff + PART_MASS] = partSys.s1[pOff + PART_MASS_VEL];
    partSys.s2dot[pOff + PART_DIAM] = 0;
    partSys.s2dot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.s2dot[pOff + PART_AGE] = .08;
    //		partSys.sdot[pOff + PART_AGE]      = (i/partSys.partCount) * 10;  // # of frame-times since creation/initialization
    partSys.s2dot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.s2dot[pOff + PART_MASS_VEL] = partSys.s1[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
    partSys.s2dot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.s2dot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.s2dot[pOff + PART_G_VEL] = partSys.s1[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
    partSys.s2dot[pOff + PART_B_VEL] = partSys.s1[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
    partSys.s2dot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.s2dot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.s2dot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu

//				console.log('Z acc = ');
//				console.log(partSys.s0[pOff + PART_ZVEL]);
}

VBObox5.prototype.veloComputer = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s1[pOff + PART_XVEL] = partSys.s0[pOff + PART_XVEL] + (partSys.s2dot[pOff + PART_XVEL] + partSys.sdot[pOff + PART_XVEL]) / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_YVEL] = partSys.s0[pOff + PART_YVEL] + (partSys.s2dot[pOff + PART_YVEL] + partSys.sdot[pOff + PART_YVEL]) / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_ZVEL] = partSys.s0[pOff + PART_ZVEL] + (partSys.s2dot[pOff + PART_ZVEL] + partSys.sdot[pOff + PART_ZVEL]) / ((timeStep * timeStep) * 2);

    partSys.s1[pOff + PART_R_VEL] = partSys.s0[pOff + PART_R_VEL] + (partSys.s2dot[pOff + PART_R_VEL] + partSys.sdot[pOff + PART_R_VEL]) / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_G_VEL] = partSys.s0[pOff + PART_G_VEL] + (partSys.s2dot[pOff + PART_G_VEL] + partSys.sdot[pOff + PART_G_VEL]) / ((timeStep * timeStep) * 2);
    partSys.s1[pOff + PART_B_VEL] = partSys.s0[pOff + PART_B_VEL] + (partSys.s2dot[pOff + PART_B_VEL] + partSys.sdot[pOff + PART_B_VEL]) / ((timeStep * timeStep) * 2);

    //partSys.s1[pOff + PART_MASS_VEL] = (partSys.s2dot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS];
    //partSys.s1[pOff + PART_AGE] = partSys.s2dot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];
}


/*
VBObox5.prototype.halfdotFinder = function(mygl, i, partSys){
				var pOff = i*PART_MAXVAR;
				partSys.sdot[pOff + PART_XPOS] = partSys.s0[pOff + PART_XVEL]/2;		// 0.0 <= randomRound() < 1.0
				partSys.sdot[pOff + PART_YPOS] = partSys.s0[pOff + PART_YVEL]/2;
				partSys.sdot[pOff + PART_ZPOS] = partSys.s0[pOff + PART_ZVEL]/2;
	//			console.log(partSys.s0[pOff + PART_X_FTOT]);
				partSys.sdot[pOff + PART_XVEL] = partSys.s0[pOff + PART_X_FTOT] * partSys.s1[pOff + PART_MASS]/2;
				partSys.sdot[pOff + PART_YVEL] = partSys.s0[pOff + PART_Y_FTOT] * partSys.s1[pOff + PART_MASS]/2;
				partSys.sdot[pOff + PART_ZVEL] = partSys.s0[pOff + PART_Z_FTOT] * partSys.s1[pOff + PART_MASS]/2;
				partSys.sdot[pOff + PART_X_FTOT] = 0.0;
				partSys.sdot[pOff + PART_Y_FTOT] = 0.0;
				partSys.sdot[pOff + PART_Z_FTOT] = 0.0;
				partSys.sdot[pOff + PART_R] = 0;
				partSys.sdot[pOff + PART_G] = partSys.s0[pOff + PART_G_VEL]/2;
				partSys.sdot[pOff + PART_B] = partSys.s0[pOff + PART_B_VEL]/2;
				partSys.sdot[pOff + PART_MASS] = partSys.s0[pOff + PART_MASS_VEL]/2;
				partSys.sdot[pOff + PART_DIAM] = 0;
				partSys.sdot[pOff + PART_RENDMODE] =0; // 0,1,2 or 3
				partSys.sdot[pOff + PART_AGE] = .08/2;
				partSys.sdot[pOff + PART_CHARGE]   =0;  // for electrostatic repulsion/attraction
				partSys.sdot[pOff + PART_MASS_VEL] = partSys.s0[pOff + PART_MASS_FTOT]/2;  // time-rate-of-change of mass.
				partSys.sdot[pOff + PART_MASS_FTOT]=0;  // force-accumulator for mass-change
				partSys.sdot[pOff + PART_R_VEL]    =0;  // time-rate-of-change of color:red
				partSys.sdot[pOff + PART_G_VEL]    = partSys.s0[pOff + PART_G_FTOT]/2;  // time-rate-of-change of color:grn
				partSys.sdot[pOff + PART_B_VEL]    = partSys.s0[pOff + PART_B_FTOT]/2; // time-rate-of-change of color:blu
				partSys.sdot[pOff + PART_R_FTOT]   =0;  // force-accumulator for color-change: red
				partSys.sdot[pOff + PART_G_FTOT]   =0;  // force-accumulator for color-change: grn
				partSys.sdot[pOff + PART_B_FTOT]   =0;  // force-accumulator for color-change: blu

}*/
/*
VBObox5.prototype.midPointFinder = function(mygl, i, partSys){
				var pOff = i*PART_MAXVAR;
				partSys.sM[pOff + PART_XPOS] = (partSys.sdot[pOff + PART_XPOS] / timeStep) + partSys.s0[pOff + PART_XPOS];
				partSys.sM[pOff + PART_YPOS] = (partSys.sdot[pOff + PART_YPOS] / timeStep) + partSys.s0[pOff + PART_YPOS];
				partSys.sM[pOff + PART_ZPOS] = (partSys.sdot[pOff + PART_ZPOS] / timeStep) + partSys.s0[pOff + PART_ZPOS];
				partSys.sM[pOff + PART_XVEL] = (partSys.sdot[pOff + PART_XVEL] / timeStep) + partSys.s0[pOff + PART_XVEL];
				partSys.sM[pOff + PART_YVEL] = (partSys.sdot[pOff + PART_YVEL] / timeStep) + partSys.s0[pOff + PART_YVEL];
				partSys.sM[pOff + PART_ZVEL] = (partSys.sdot[pOff + PART_ZVEL] / timeStep) + partSys.s0[pOff + PART_ZVEL];
				partSys.sM[pOff + PART_R] = (partSys.sdot[pOff + PART_R] / timeStep) + partSys.s0[pOff + PART_R];
				partSys.sM[pOff + PART_G] = (partSys.sdot[pOff + PART_G] / timeStep) + partSys.s0[pOff + PART_G];
				partSys.sM[pOff + PART_B] = (partSys.sdot[pOff + PART_B] / timeStep) + partSys.s0[pOff + PART_B];
				partSys.sM[pOff + PART_R_VEL] = (partSys.sdot[pOff + PART_R_VEL] / timeStep) + partSys.s0[pOff + PART_R_VEL];
				partSys.sM[pOff + PART_G_VEL] = (partSys.sdot[pOff + PART_G_VEL] / timeStep) + partSys.s0[pOff + PART_G_VEL];
				partSys.sM[pOff + PART_B_VEL] = (partSys.sdot[pOff + PART_B_VEL] / timeStep) + partSys.s0[pOff + PART_B_VEL];
				partSys.sM[pOff + PART_MASS] = (partSys.sdot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS];
				partSys.sM[pOff + PART_AGE] = partSys.sdot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];
				
}

VBObox5.prototype.applyMidForces = function(mygl, i, partSys) {

			var rad = 0;
			
				var pOff = i*PART_MAXVAR;
				partSys.sM[pOff + PART_X_FTOT] = 0;
				partSys.sM[pOff + PART_Y_FTOT] = 0;
				partSys.sM[pOff + PART_Z_FTOT] = 0;
				partSys.cforcer[F_MOUSE](pOff, partSys.sM);
				partSys.cforcer[F_GRAV_E](pOff, partSys.sM);
				partSys.cforcer[F_GRAV_P](pOff, partSys.sM);
				partSys.cforcer[F_WIND](pOff, partSys.sM);
				partSys.cforcer[F_BUBBLE](pOff, partSys.sM);
				partSys.cforcer[F_DRAG](pOff, partSys.sM);
				
				partSys.cforcer[F_SPRING](pOff, partSys.sM);
				partSys.cforcer[F_SPH](pOff, partSys.sM);
				partSys.cforcer[F_CHARGE](pOff, partSys.sM);
				partSys.cforcer[F_FIRE](pOff, partSys.sM);
			*/
/*const F_NONE      0       // Non-existent force: ignore this CForcer object
const F_MOUSE     1       // Spring-like connection to the mouse cursor; lets
            // you 'grab' and 'wiggle' one particle(or several).
const F_GRAV_E    2       // Earth-gravity: pulls all particles 'downward'.
const F_GRAV_P    3       // Planetary-gravity; particle-pair (e0,e1) attract
            // each other with force== grav* mass0*mass1/ dist^2
const F_WIND      4       // Blowing-wind-like force-field;fcn of 3D position
const F_BUBBLE    5       // Constant inward force towards centerpoint if
            // particle is > max_radius away from centerpoint.
const F_DRAG      6       // Viscous drag -- proportional to neg. velocity.
const F_SPRING    7       // ties together 2 particles; distance sets force
const F_SPH 8       // a big collection of identical springs; lets you
            // make cloth & rubbery shapes as one force-making
            // object, instead of many many F_SPRING objects.
const F_CHARGE    9       // attract/repel by charge and inverse distance;
            // applies to all charged particles.
const F_MAXKINDS  10      // 'max' is always the LAST name in our list;
*/
/*
}

VBObox5.prototype.midDotFinder = function(mygl, i, partSys){
				var pOff = i*PART_MAXVAR;
				partSys.s2dot[pOff + PART_XPOS] = partSys.sM[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
				partSys.s2dot[pOff + PART_YPOS] = partSys.sM[pOff + PART_YVEL];
				partSys.s2dot[pOff + PART_ZPOS] = partSys.sM[pOff + PART_ZVEL];
				partSys.s2dot[pOff + PART_XVEL] = partSys.sM[pOff + PART_X_FTOT] * partSys.sM[pOff + PART_MASS];
				partSys.s2dot[pOff + PART_YVEL] = partSys.sM[pOff + PART_Y_FTOT] * partSys.sM[pOff + PART_MASS];
				partSys.s2dot[pOff + PART_ZVEL] = partSys.sM[pOff + PART_Z_FTOT] * partSys.sM[pOff + PART_MASS];
				partSys.s2dot[pOff + PART_X_FTOT] = 0.0;
				partSys.s2dot[pOff + PART_Y_FTOT] = 0.0;
				partSys.s2dot[pOff + PART_Z_FTOT] = 0.0;
				partSys.s2dot[pOff + PART_R] = 0;
				partSys.s2dot[pOff + PART_G] = partSys.sM[pOff + PART_G_VEL];
				partSys.s2dot[pOff + PART_B] = partSys.sM[pOff + PART_B_VEL];
				partSys.s2dot[pOff + PART_MASS] = partSys.sM[pOff + PART_MASS_VEL];
				partSys.s2dot[pOff + PART_DIAM] = 0;
				partSys.s2dot[pOff + PART_RENDMODE] =0; // 0,1,2 or 3
				partSys.s2dot[pOff + PART_AGE] = .08;
		//		partSys.s2dot[pOff + PART_AGE]      = (i/partSys.partCount) * 10;  // # of frame-times since creation/initialization
				partSys.s2dot[pOff + PART_CHARGE]   =0;  // for electrostatic repulsion/attraction
				partSys.s2dot[pOff + PART_MASS_VEL] = partSys.sM[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
				partSys.s2dot[pOff + PART_MASS_FTOT]=0;  // force-accumulator for mass-change
				partSys.s2dot[pOff + PART_R_VEL]    =0;  // time-rate-of-change of color:red
				partSys.s2dot[pOff + PART_G_VEL]    = partSys.sM[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
				partSys.s2dot[pOff + PART_B_VEL]    = partSys.sM[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
				partSys.s2dot[pOff + PART_R_FTOT]   =0;  // force-accumulator for color-change: red
				partSys.s2dot[pOff + PART_G_FTOT]   =0;  // force-accumulator for color-change: grn
				partSys.s2dot[pOff + PART_B_FTOT]   =0;  // force-accumulator for color-change: blu

}

VBObox5.prototype.render = function(mygl, i, partSys){
				var pOff = i*PART_MAXVAR;
				partSys.s1[pOff + PART_XPOS] = (partSys.s2dot[pOff + PART_XPOS] / timeStep) + partSys.s0[pOff + PART_XPOS];
				partSys.s1[pOff + PART_YPOS] = (partSys.s2dot[pOff + PART_YPOS] / timeStep) + partSys.s0[pOff + PART_YPOS];
				partSys.s1[pOff + PART_ZPOS] = (partSys.s2dot[pOff + PART_ZPOS] / timeStep) + partSys.s0[pOff + PART_ZPOS];
				partSys.s1[pOff + PART_XVEL] = (partSys.s2dot[pOff + PART_XVEL] / timeStep) + partSys.s0[pOff + PART_XVEL];
				partSys.s1[pOff + PART_YVEL] = (partSys.s2dot[pOff + PART_YVEL] / timeStep) + partSys.s0[pOff + PART_YVEL];
				partSys.s1[pOff + PART_ZVEL] = (partSys.s2dot[pOff + PART_ZVEL] / timeStep) + partSys.s0[pOff + PART_ZVEL];
				partSys.s1[pOff + PART_R] = (partSys.s2dot[pOff + PART_R] /timeStep) + partSys.s0[pOff + PART_R];
				partSys.s1[pOff + PART_G] = (partSys.s2dot[pOff + PART_G] / timeStep) + partSys.s0[pOff + PART_G];
				partSys.s1[pOff + PART_B] = (partSys.s2dot[pOff + PART_B] / timeStep) + partSys.s0[pOff + PART_B];
				partSys.s1[pOff + PART_R_VEL] = (partSys.s2dot[pOff + PART_R_VEL] / timeStep) + partSys.s0[pOff + PART_R_VEL];
				partSys.s1[pOff + PART_G_VEL] = (partSys.s2dot[pOff + PART_G_VEL] / timeStep) + partSys.s0[pOff + PART_G_VEL];
				partSys.s1[pOff + PART_B_VEL] = (partSys.s2dot[pOff + PART_B_VEL] / timeStep) + partSys.s0[pOff + PART_B_VEL];
				partSys.s1[pOff + PART_MASS] = (partSys.s2dot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS];
				partSys.s1[pOff + PART_AGE] = partSys.s2dot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];
				
	//			if(pOff == 0) console.log(partSys.s1[pOff + PART_XPOS]);
		//		console.log(partSys.s1[pOff + PART_ZVEL]);
}
*/
VBObox5.prototype.PartSys_constrain1 = function (i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.c0[0](i, partSys.s1, partSys.s0);
    partSys.c0[1](pOff, partSys.s1);
}

VBObox5.prototype.PartSys_render = function (myGL, partSys) {

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
    // 	myGL.uniform1i(this.u_runModeID, myRunMode);		// run/step/pause the particle system

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
    myGL.bufferSubData(myGL.ARRAY_BUFFER, 0, partSys.s0);
    myGL.drawArrays(myGL.POINTS, 0, partSys.partCount);

}


VBObox5.prototype.reset = function (partSys) {
    /*for (i = 0; i < partSys.partCount; i++){
                if(s0[i* PART_MAXVAR + PART_XVEL] > 0.0) s0[i* PART_MAXVAR + PART_XVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_XVEL] -= INIT_VEL;
                if(s0[i* PART_MAXVAR + PART_YVEL] > 0.0) s0[i* PART_MAXVAR + PART_YVEL] += INIT_VEL*.8; else s0[i* PART_MAXVAR + PART_YVEL] -= INIT_VEL*.8;
                if(s0[i* PART_MAXVAR + PART_ZVEL] > 0.0) s0[i* PART_MAXVAR + PART_ZVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_ZVEL] -= INIT_VEL;
                }*/
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s0[pOff + PART_XVEL] > 0) {
            partSys.s0[pOff + PART_XVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_XVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_YVEL] > 0) {
            partSys.s0[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_ZVEL] > 0) {
            partSys.s0[pOff + PART_ZVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_ZVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
    }
}

VBObox5.prototype.resetWind = function (partSys) {
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s0[pOff + PART_XVEL] > 0) {
            partSys.s0[pOff + PART_XVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_XVEL] -= (0.3 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_YVEL] > 0) {
            partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
    }
}

VBObox5.prototype.resetSpring = function (partSys) {
    for (var i = 1; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        partSys.s0[pOff + PART_ZPOS] -= springRun * (.002 * Math.pow(i, 2));
        /*			}
                else partSys.s0[pOff + PART_XVEL] -= (0.3 + 0.8*Math.random())*INIT_VEL;
                if(  partSys.s0[pOff + PART_YVEL] > 0) {
                       partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8*Math.random())*INIT_VEL;
                    }
                else partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8*Math.random())*INIT_VEL;
            }*/
    }
}


function VBObox6() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox' object  that holds all data and fcns 
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate set of shaders.
    this.VERT_SRC =	//--------------------- VERTEX SHADER source code
        'precision mediump float;\n' +				// req'd in OpenGL ES if we use 'float'
        //
        'uniform   int u_runMode1; \n' +					// particle system state:
        // 0=reset; 1= pause; 2=step; 3=run
        // 'uniform	 vec4 u_ballShift1; \n' +			// single bouncy-ball's movement
        'uniform   mat4 u_MvpMat1; \n' +    //camera position
        'attribute vec3 a_Pos1;\n' +
        'attribute vec3 a_Colr1;\n' +
        'varying   vec4 v_Colr1; \n' +
        'attribute float a_diam1; \n' +
        'void main() {\n' +
        '  gl_PointSize = a_diam1;\n' +
        '	 gl_Position = u_MvpMat1 * (vec4(a_Pos1, 1)); \n' +// + u_ballShift1)

        // Let u_runMode determine particle color:
        //  '  if(u_runMode1 == 0) { \n' +
        //	'	   v_Colr1 = vec4(1.0, 0.0, 0.0, 1.0);	\n' +		// red: 0==reset
        //	'  	 } \n' +
        //	'  else if(u_runMode1 == 1) {  \n' +
        //	'    v_Colr1 = vec4(1.0, 1.0, 0.0, 1.0); \n' +	// yellow: 1==pause
        //	'    }  \n' +
        //	'  else if(u_runMode1 == 2) { \n' +
        //	'    v_Colr1 = vec4(1.0, 1.0, 1.0, 1.0); \n' +	// white: 2==step
        //'    } \n' +
        //'  else { \n' +
        'v_Colr1 = vec4(a_Colr1, 1.0); \n' +	// green: >3==run
        //	'		 } \n' +
        '} \n';

    this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code
        'precision mediump float;\n' +
        'varying vec4 v_Colr1; \n' +
        'void main() {\n' +
        '  float dist = distance(gl_PointCoord, vec2(0.5, 0.5)); \n' +
        //  '  if(dist < 0.5) { \n' +
        '  	gl_FragColor = vec4((1.0-2.0*dist)*v_Colr1.rgb, 1.0);\n' +
        //	'  } else { discard; }\n' +
        '}\n';


    this.vboLoc;										// Vertex Buffer Object location# on the GPU
    this.vertices;
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


VBObox6.prototype.init = function (myGL) {
//=============================================================================
// Create, compile, link this VBObox object's shaders to an executable 'program'
// ready for use in the GPU.  Create and fill a Float32Array that holds all VBO 
// vertices' values; create a new VBO on the GPU and fill it with those values. 
// Find the GPU location of	all our shaders' attribute- and uniform-variables; 
// assign the correct portions of VBO contents as the data source for each 
// attribute, and transfer current values to the GPU for each uniform variable.
// (usually called only once, within main()) 

//partSys.s1 = new Float32Array(partSys.partCount*PART_MAXVAR);
    this.FSIZE = partSys1.s0.BYTES_PER_ELEMENT;	// memory needed to store an s0 array element.

// this.PartSys_init();
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
        partSys1.s1, 		// JavaScript Float32Array
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


}


VBObox6.prototype.adjust1 = function (mygl, partSys) {
    mygl.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    for (var i = 0; i < partSys.partCount; i++) {
        this.applyForces(mygl, i, partSys);

        this.dotFinder(mygl, i, partSys);
        if (partSys.first == 0)
            this.oldPointFinder(mygl, i, partSys);
        this.apply2Forces(mygl, i, partSys);
        this.Dot2Finder(mygl, i, partSys);
        this.render(mygl, i, partSys);
        this.PartSys_constrain1(i, partSys);
    }
    partSys.buffer = partSys.s1.slice();
    partSys.s1dot2 = partSys.s2.slice();
    partSys.s0 = partSys.s1.slice();
    partSys.s1 = partSys.s1dot2.slice();
    partSys.s2 = partSys.buffer.slice();
    this.vertices = partSys.s1dot2;
    this.PartSys_render(gl, partSys);
}


VBObox6.prototype.applyForces = function (mygl, i, partSys) {

    var rad = 0;
    var pOff = i * PART_MAXVAR;
    partSys.s1[pOff + PART_X_FTOT] = 0;
    partSys.s1[pOff + PART_Y_FTOT] = 0;
    partSys.s1[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.s1);
    partSys.cforcer[F_GRAV_E](pOff, partSys.s1);
    partSys.cforcer[F_GRAV_P](pOff, partSys.s1);
    partSys.cforcer[F_WIND](pOff, partSys.s1);
    partSys.cforcer[F_BUBBLE](pOff, partSys.s1);
    partSys.cforcer[F_DRAG](pOff, partSys.s1);
    partSys.cforcer[F_SPRING](pOff, partSys.s1);
    partSys.cforcer[F_SPH](pOff, partSys.s1);
    partSys.cforcer[F_CHARGE](pOff, partSys.s1);
    partSys.cforcer[F_FIRE](pOff, partSys.s1);

    /*const F_NONE      0       // Non-existent force: ignore this CForcer object
const F_MOUSE     1       // Spring-like connection to the mouse cursor; lets
                // you 'grab' and 'wiggle' one particle(or several).
const F_GRAV_E    2       // Earth-gravity: pulls all particles 'downward'.
const F_GRAV_P    3       // Planetary-gravity; particle-pair (e0,e1) attract
                // each other with force== grav* mass0*mass1/ dist^2
const F_WIND      4       // Blowing-wind-like force-field;fcn of 3D position
const F_BUBBLE    5       // Constant inward force towards centerpoint if
                // particle is > max_radius away from centerpoint.
const F_DRAG      6       // Viscous drag -- proportional to neg. velocity.
const F_SPRING    7       // ties together 2 particles; distance sets force
const F_SPH 8       // a big collection of identical springs; lets you
                // make cloth & rubbery shapes as one force-making
                // object, instead of many many F_SPRING objects.
const F_CHARGE    9       // attract/repel by charge and inverse distance;
                // applies to all charged particles.
const F_MAXKINDS  10      // 'max' is always the LAST name in our list;
*/

}

VBObox6.prototype.dotFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.sdot[pOff + PART_XPOS] = partSys.s1[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
    partSys.sdot[pOff + PART_YPOS] = partSys.s1[pOff + PART_YVEL];
    partSys.sdot[pOff + PART_ZPOS] = partSys.s1[pOff + PART_ZVEL];
    partSys.sdot[pOff + PART_XVEL] = partSys.s1[pOff + PART_X_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.sdot[pOff + PART_YVEL] = partSys.s1[pOff + PART_Y_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.sdot[pOff + PART_ZVEL] = partSys.s1[pOff + PART_Z_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.sdot[pOff + PART_X_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Y_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Z_FTOT] = 0.0;
    partSys.sdot[pOff + PART_R] = 0;
    partSys.sdot[pOff + PART_G] = partSys.s1[pOff + PART_G_VEL];
    partSys.sdot[pOff + PART_B] = partSys.s1[pOff + PART_B_VEL];
    partSys.sdot[pOff + PART_MASS] = partSys.s1[pOff + PART_MASS_VEL];
    partSys.sdot[pOff + PART_DIAM] = 0;
    partSys.sdot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.sdot[pOff + PART_AGE] = .08;
    //		partSys.sdot[pOff + PART_AGE]      = (i/partSys.partCount) * 10;  // # of frame-times since creation/initialization
    partSys.sdot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.sdot[pOff + PART_MASS_VEL] = partSys.s1[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
    partSys.sdot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.sdot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.sdot[pOff + PART_G_VEL] = partSys.s1[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
    partSys.sdot[pOff + PART_B_VEL] = partSys.s1[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
    partSys.sdot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.sdot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.sdot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu

//				console.log('Z acc = ');
//				console.log(partSys.s0[pOff + PART_ZVEL]);
}

VBObox6.prototype.oldPointFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s0[pOff + PART_XPOS] = partSys.s1[pOff + PART_XPOS] - (partSys.sdot[pOff + PART_XPOS] / timeStep);
    partSys.s0[pOff + PART_YPOS] = partSys.s1[pOff + PART_YPOS] - (partSys.sdot[pOff + PART_YPOS] / timeStep);
    partSys.s0[pOff + PART_ZPOS] = partSys.s1[pOff + PART_ZPOS] - (partSys.sdot[pOff + PART_ZPOS] / timeStep);
    partSys.s0[pOff + PART_XVEL] = partSys.s1[pOff + PART_XVEL] - (partSys.sdot[pOff + PART_XVEL] / timeStep);
    partSys.s0[pOff + PART_YVEL] = partSys.s1[pOff + PART_YVEL] - (partSys.sdot[pOff + PART_YVEL] / timeStep);
    partSys.s0[pOff + PART_ZVEL] = partSys.s1[pOff + PART_ZVEL] - (partSys.sdot[pOff + PART_ZVEL] / timeStep);
    partSys.s0[pOff + PART_R] = partSys.s1[pOff + PART_R] - (partSys.sdot[pOff + PART_R] / timeStep);
    partSys.s0[pOff + PART_G] = partSys.s1[pOff + PART_G] - (partSys.sdot[pOff + PART_G] / timeStep);
    partSys.s0[pOff + PART_B] = partSys.s1[pOff + PART_B] - (partSys.sdot[pOff + PART_B] / timeStep);
    partSys.s0[pOff + PART_R_VEL] = partSys.s1[pOff + PART_R_VEL] - (partSys.sdot[pOff + PART_R_VEL] / timeStep);
    partSys.s0[pOff + PART_G_VEL] = partSys.s1[pOff + PART_G_VEL] - (partSys.sdot[pOff + PART_G_VEL] / timeStep);
    partSys.s0[pOff + PART_B_VEL] = partSys.s1[pOff + PART_B_VEL] - (partSys.sdot[pOff + PART_B_VEL] / timeStep);
    partSys.s0[pOff + PART_MASS] = partSys.s1[pOff + PART_MASS] - (partSys.sdot[pOff + PART_MASS] / timeStep);
    partSys.s0[pOff + PART_AGE] = partSys.s1[pOff + PART_AGE] - partSys.sdot[pOff + PART_AGE];
    partSys.first = 1;
}

VBObox6.prototype.apply2Forces = function (mygl, i, partSys) {

    var rad = 0;
    var pOff = i * PART_MAXVAR;
    partSys.s0[pOff + PART_X_FTOT] = 0;
    partSys.s0[pOff + PART_Y_FTOT] = 0;
    partSys.s0[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_E](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_P](pOff, partSys.s0);
    partSys.cforcer[F_WIND](pOff, partSys.s0);
    partSys.cforcer[F_BUBBLE](pOff, partSys.s0);
    partSys.cforcer[F_DRAG](pOff, partSys.s0);
    partSys.cforcer[F_SPRING](pOff, partSys.s0);
    partSys.cforcer[F_SPH](pOff, partSys.s0);
    partSys.cforcer[F_CHARGE](pOff, partSys.s0);
    partSys.cforcer[F_FIRE](pOff, partSys.s0);

    /*const F_NONE      0       // Non-existent force: ignore this CForcer object
const F_MOUSE     1       // Spring-like connection to the mouse cursor; lets
                // you 'grab' and 'wiggle' one particle(or several).
const F_GRAV_E    2       // Earth-gravity: pulls all particles 'downward'.
const F_GRAV_P    3       // Planetary-gravity; particle-pair (e0,e1) attract
                // each other with force== grav* mass0*mass1/ dist^2
const F_WIND      4       // Blowing-wind-like force-field;fcn of 3D position
const F_BUBBLE    5       // Constant inward force towards centerpoint if
                // particle is > max_radius away from centerpoint.
const F_DRAG      6       // Viscous drag -- proportional to neg. velocity.
const F_SPRING    7       // ties together 2 particles; distance sets force
const F_SPH 8       // a big collection of identical springs; lets you
                // make cloth & rubbery shapes as one force-making
                // object, instead of many many F_SPRING objects.
const F_CHARGE    9       // attract/repel by charge and inverse distance;
                // applies to all charged particles.
const F_MAXKINDS  10      // 'max' is always the LAST name in our list;
*/
}

VBObox6.prototype.Dot2Finder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s2dot[pOff + PART_XPOS] = partSys.s0[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
    partSys.s2dot[pOff + PART_YPOS] = partSys.s0[pOff + PART_YVEL];
    partSys.s2dot[pOff + PART_ZPOS] = partSys.s0[pOff + PART_ZVEL];

    partSys.s2dot[pOff + PART_XVEL] = partSys.s0[pOff + PART_X_FTOT] * partSys.s0[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_YVEL] = partSys.s0[pOff + PART_Y_FTOT] * partSys.s0[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_ZVEL] = partSys.s0[pOff + PART_Z_FTOT] * partSys.s0[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_X_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_Y_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_Z_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_R] = 0;
    partSys.s2dot[pOff + PART_G] = partSys.s0[pOff + PART_G_VEL];
    partSys.s2dot[pOff + PART_B] = partSys.s0[pOff + PART_B_VEL];
    partSys.s2dot[pOff + PART_MASS] = partSys.s0[pOff + PART_MASS_VEL];
    partSys.s2dot[pOff + PART_DIAM] = 0;
    partSys.s2dot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.s2dot[pOff + PART_AGE] = .08;
    //		partSys.sdot[pOff + PART_AGE]      = (i/partSys.partCount) * 10;  // # of frame-times since creation/initialization
    partSys.s2dot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.s2dot[pOff + PART_MASS_VEL] = partSys.s0[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
    partSys.s2dot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.s2dot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.s2dot[pOff + PART_G_VEL] = partSys.s0[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
    partSys.s2dot[pOff + PART_B_VEL] = partSys.s0[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
    partSys.s2dot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.s2dot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.s2dot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
}


VBObox6.prototype.render = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s2[pOff + PART_XPOS] = (partSys.sdot[pOff + PART_XPOS] / (timeStep * (2 / 3))) + partSys.s1[pOff + PART_XPOS] - (partSys.s2dot[pOff + PART_XPOS] / (timeStep * 2));
    partSys.s2[pOff + PART_YPOS] = (partSys.sdot[pOff + PART_YPOS] / (timeStep * (2 / 3))) + partSys.s1[pOff + PART_YPOS] - (partSys.s2dot[pOff + PART_YPOS] / (timeStep * 2));
    partSys.s2[pOff + PART_ZPOS] = (partSys.sdot[pOff + PART_ZPOS] / (timeStep * (2 / 3))) + partSys.s1[pOff + PART_ZPOS] - (partSys.s2dot[pOff + PART_ZPOS] / (timeStep * 2));
    partSys.s2[pOff + PART_XVEL] = (partSys.sdot[pOff + PART_XVEL] / (timeStep * (2 / 3))) + partSys.s1[pOff + PART_XVEL] - (partSys.s2dot[pOff + PART_XVEL] / (timeStep * 2));
    partSys.s2[pOff + PART_YVEL] = (partSys.sdot[pOff + PART_YVEL] / (timeStep * (2 / 3))) + partSys.s1[pOff + PART_YVEL] - (partSys.s2dot[pOff + PART_YVEL] / (timeStep * 2));
    partSys.s2[pOff + PART_ZVEL] = (partSys.sdot[pOff + PART_ZVEL] / (timeStep * (2 / 3))) + partSys.s1[pOff + PART_ZVEL] - (partSys.s2dot[pOff + PART_ZVEL] / (timeStep * 2));
    partSys.s2[pOff + PART_R] = (partSys.sdot[pOff + PART_R] / (timeStep * (2 / 3))) + partSys.s0[pOff + PART_R] - (partSys.s2dot[pOff + PART_R] / (timeStep * 2));
    partSys.s2[pOff + PART_G] = (partSys.sdot[pOff + PART_G] / (timeStep * (2 / 3))) + partSys.s0[pOff + PART_G] - (partSys.s2dot[pOff + PART_G] / (timeStep * 2));
    partSys.s2[pOff + PART_B] = (partSys.sdot[pOff + PART_B] / (timeStep * (2 / 3))) + partSys.s0[pOff + PART_B] - (partSys.s2dot[pOff + PART_B] / (timeStep * 2));
    partSys.s2[pOff + PART_R_VEL] = (partSys.sdot[pOff + PART_R_VEL] / (timeStep * (2 / 3))) + partSys.s0[pOff + PART_R_VEL] - (partSys.s2dot[pOff + PART_R_VEL] / (timeStep * 2));
    partSys.s2[pOff + PART_G_VEL] = (partSys.sdot[pOff + PART_G_VEL] / (timeStep * (2 / 3))) + partSys.s0[pOff + PART_G_VEL] - (partSys.s2dot[pOff + PART_G_VEL] / (timeStep * 2));
    partSys.s2[pOff + PART_B_VEL] = (partSys.sdot[pOff + PART_B_VEL] / (timeStep * (2 / 3))) + partSys.s0[pOff + PART_B_VEL] - (partSys.s2dot[pOff + PART_B_VEL] / (timeStep * 2));
    partSys.s2[pOff + PART_MASS] = (partSys.sdot[pOff + PART_MASS] / (timeStep * (2 / 3))) + partSys.s0[pOff + PART_MASS] - (partSys.s2dot[pOff + PART_MASS] / (timeStep * 2));
    partSys.s2[pOff + PART_AGE] = partSys.sdot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];

    //			if(pOff == 0) console.log(partSys.s1[pOff + PART_XPOS]);
    //		console.log(partSys.s1[pOff + PART_ZVEL]);
}

VBObox6.prototype.PartSys_constrain1 = function (i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.c0[0](i, partSys.s2, partSys.s1);
    partSys.c0[1](pOff, partSys.s2);
}

VBObox6.prototype.PartSys_render = function (myGL, partSys) {

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
    // 	myGL.uniform1i(this.u_runModeID, myRunMode);		// run/step/pause the particle system

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
    myGL.bufferSubData(myGL.ARRAY_BUFFER, 0, partSys.s1);
    myGL.drawArrays(myGL.POINTS, 0, partSys.partCount);

}


VBObox6.prototype.reset = function (partSys) {
    /*for (i = 0; i < partSys.partCount; i++){
                if(s0[i* PART_MAXVAR + PART_XVEL] > 0.0) s0[i* PART_MAXVAR + PART_XVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_XVEL] -= INIT_VEL;
                if(s0[i* PART_MAXVAR + PART_YVEL] > 0.0) s0[i* PART_MAXVAR + PART_YVEL] += INIT_VEL*.8; else s0[i* PART_MAXVAR + PART_YVEL] -= INIT_VEL*.8;
                if(s0[i* PART_MAXVAR + PART_ZVEL] > 0.0) s0[i* PART_MAXVAR + PART_ZVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_ZVEL] -= INIT_VEL;
                }*/
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s1[pOff + PART_XVEL] > 0) {
            partSys.s1[pOff + PART_XVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s1[pOff + PART_XVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s1[pOff + PART_YVEL] > 0) {
            partSys.s1[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s1[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s1[pOff + PART_ZVEL] > 0) {
            partSys.s1[pOff + PART_ZVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s1[pOff + PART_ZVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
    }
}

VBObox6.prototype.resetWind = function (partSys) {
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s1[pOff + PART_XVEL] > 0) {
            partSys.s1[pOff + PART_XVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s1[pOff + PART_XVEL] -= (0.3 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s1[pOff + PART_YVEL] > 0) {
            partSys.s1[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s1[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
    }
}

VBObox6.prototype.resetSpring = function (partSys) {
    for (var i = 1; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        partSys.s1[pOff + PART_ZPOS] -= springRun * (.002 * Math.pow(i, 2));
        /*			}
                else partSys.s0[pOff + PART_XVEL] -= (0.3 + 0.8*Math.random())*INIT_VEL;
                if(  partSys.s0[pOff + PART_YVEL] > 0) {
                       partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8*Math.random())*INIT_VEL;
                    }
                else partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8*Math.random())*INIT_VEL;
            }*/
    }
}

function VBObox7() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox' object  that holds all data and fcns 
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate set of shaders.
    this.VERT_SRC =	//--------------------- VERTEX SHADER source code
        'precision mediump float;\n' +				// req'd in OpenGL ES if we use 'float'
        //
        'uniform   int u_runMode1; \n' +					// particle system state:
        // 0=reset; 1= pause; 2=step; 3=run
        // 'uniform	 vec4 u_ballShift1; \n' +			// single bouncy-ball's movement
        'uniform   mat4 u_MvpMat1; \n' +    //camera position
        'attribute vec3 a_Pos1;\n' +
        'attribute vec3 a_Colr1;\n' +
        'varying   vec4 v_Colr1; \n' +
        'attribute float a_diam1; \n' +
        'void main() {\n' +
        '  gl_PointSize = a_diam1;\n' +
        '	 gl_Position = u_MvpMat1 * (vec4(a_Pos1, 1)); \n' +// + u_ballShift1)

        // Let u_runMode determine particle color:
        //  '  if(u_runMode1 == 0) { \n' +
        //	'	   v_Colr1 = vec4(1.0, 0.0, 0.0, 1.0);	\n' +		// red: 0==reset
        //	'  	 } \n' +
        //	'  else if(u_runMode1 == 1) {  \n' +
        //	'    v_Colr1 = vec4(1.0, 1.0, 0.0, 1.0); \n' +	// yellow: 1==pause
        //	'    }  \n' +
        //	'  else if(u_runMode1 == 2) { \n' +
        //	'    v_Colr1 = vec4(1.0, 1.0, 1.0, 1.0); \n' +	// white: 2==step
        //'    } \n' +
        //'  else { \n' +
        'v_Colr1 = vec4(a_Colr1, 1.0); \n' +	// green: >3==run
        //	'		 } \n' +
        '} \n';

    this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code
        'precision mediump float;\n' +
        'varying vec4 v_Colr1; \n' +
        'void main() {\n' +
        '  float dist = distance(gl_PointCoord, vec2(0.5, 0.5)); \n' +
        //  '  if(dist < 0.5) { \n' +
        '  	gl_FragColor = vec4((1.0-2.0*dist)*v_Colr1.rgb, 1.0);\n' +
        //	'  } else { discard; }\n' +
        '}\n';


    this.vboLoc;										// Vertex Buffer Object location# on the GPU
    this.vertices;
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


VBObox7.prototype.init = function (myGL) {
//=============================================================================
// Create, compile, link this VBObox object's shaders to an executable 'program'
// ready for use in the GPU.  Create and fill a Float32Array that holds all VBO 
// vertices' values; create a new VBO on the GPU and fill it with those values. 
// Find the GPU location of	all our shaders' attribute- and uniform-variables; 
// assign the correct portions of VBO contents as the data source for each 
// attribute, and transfer current values to the GPU for each uniform variable.
// (usually called only once, within main()) 

//partSys.s1 = new Float32Array(partSys.partCount*PART_MAXVAR);
    this.FSIZE = partSys1.s0.BYTES_PER_ELEMENT;	// memory needed to store an s0 array element.

// this.PartSys_init();
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
        partSys1.s0, 		// JavaScript Float32Array
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


}


VBObox7.prototype.adjust1 = function (mygl, partSys) {
    mygl.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    for (var i = 0; i < partSys.partCount; i++) {
        this.applyForces(mygl, i, partSys);
        this.halfdotFinder(mygl, i, partSys);
        this.midPointFinder(mygl, i, partSys);
        this.applyMidForces(mygl, i, partSys);
        this.midDotFinder(mygl, i, partSys);
        this.render(mygl, i, partSys);
        this.applyForces3(mygl, i, partSys);
        this.revDotFinder(mygl, i, partSys);
        this.render2(mygl, i, partSys);
        this.PartSys_constrain1(i, partSys);
    }
    partSys.buffer = partSys.s0.slice();
    partSys.s0 = partSys.s1.slice();
    partSys.s1 = partSys.buffer.slice();
    this.vertices = partSys.s0;
    this.PartSys_render(gl, partSys);

}

VBObox7.prototype.applyForces = function (mygl, i, partSys) {

    var rad = 0;
    var pOff = i * PART_MAXVAR;
    partSys.s0[pOff + PART_X_FTOT] = 0;
    partSys.s0[pOff + PART_Y_FTOT] = 0;
    partSys.s0[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_E](pOff, partSys.s0);
    partSys.cforcer[F_GRAV_P](pOff, partSys.s0);
    partSys.cforcer[F_WIND](pOff, partSys.s0);
    partSys.cforcer[F_BUBBLE](pOff, partSys.s0);
    partSys.cforcer[F_DRAG](pOff, partSys.s0);
    partSys.cforcer[F_SPRING](pOff, partSys.s0);
    partSys.cforcer[F_SPH](pOff, partSys.s0);
    partSys.cforcer[F_CHARGE](pOff, partSys.s0);
    partSys.cforcer[F_FIRE](pOff, partSys.s0);

    /*const F_NONE      0       // Non-existent force: ignore this CForcer object
const F_MOUSE     1       // Spring-like connection to the mouse cursor; lets
                // you 'grab' and 'wiggle' one particle(or several).
const F_GRAV_E    2       // Earth-gravity: pulls all particles 'downward'.
const F_GRAV_P    3       // Planetary-gravity; particle-pair (e0,e1) attract
                // each other with force== grav* mass0*mass1/ dist^2
const F_WIND      4       // Blowing-wind-like force-field;fcn of 3D position
const F_BUBBLE    5       // Constant inward force towards centerpoint if
                // particle is > max_radius away from centerpoint.
const F_DRAG      6       // Viscous drag -- proportional to neg. velocity.
const F_SPRING    7       // ties together 2 particles; distance sets force
const F_SPH 8       // a big collection of identical springs; lets you
                // make cloth & rubbery shapes as one force-making
                // object, instead of many many F_SPRING objects.
const F_CHARGE    9       // attract/repel by charge and inverse distance;
                // applies to all charged particles.
const F_MAXKINDS  10      // 'max' is always the LAST name in our list;
*/

}

VBObox7.prototype.halfdotFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.sdot[pOff + PART_XPOS] = partSys.s0[pOff + PART_XVEL] / 2;		// 0.0 <= randomRound() < 1.0
    partSys.sdot[pOff + PART_YPOS] = partSys.s0[pOff + PART_YVEL] / 2;
    partSys.sdot[pOff + PART_ZPOS] = partSys.s0[pOff + PART_ZVEL] / 2;
    //			console.log(partSys.s0[pOff + PART_X_FTOT]);
    partSys.sdot[pOff + PART_XVEL] = partSys.s0[pOff + PART_X_FTOT] * partSys.s1[pOff + PART_MASS] / 2;
    partSys.sdot[pOff + PART_YVEL] = partSys.s0[pOff + PART_Y_FTOT] * partSys.s1[pOff + PART_MASS] / 2;
    partSys.sdot[pOff + PART_ZVEL] = partSys.s0[pOff + PART_Z_FTOT] * partSys.s1[pOff + PART_MASS] / 2;
    partSys.sdot[pOff + PART_X_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Y_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Z_FTOT] = 0.0;
    partSys.sdot[pOff + PART_R] = 0;
    partSys.sdot[pOff + PART_G] = partSys.s0[pOff + PART_G_VEL] / 2;
    partSys.sdot[pOff + PART_B] = partSys.s0[pOff + PART_B_VEL] / 2;
    partSys.sdot[pOff + PART_MASS] = partSys.s0[pOff + PART_MASS_VEL] / 2;
    partSys.sdot[pOff + PART_DIAM] = 0;
    partSys.sdot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.sdot[pOff + PART_AGE] = .08 / 2;
    partSys.sdot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.sdot[pOff + PART_MASS_VEL] = partSys.s0[pOff + PART_MASS_FTOT] / 2;  // time-rate-of-change of mass.
    partSys.sdot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.sdot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.sdot[pOff + PART_G_VEL] = partSys.s0[pOff + PART_G_FTOT] / 2;  // time-rate-of-change of color:grn
    partSys.sdot[pOff + PART_B_VEL] = partSys.s0[pOff + PART_B_FTOT] / 2; // time-rate-of-change of color:blu
    partSys.sdot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.sdot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.sdot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu

}

VBObox7.prototype.midPointFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.sM[pOff + PART_XPOS] = (partSys.sdot[pOff + PART_XPOS] / timeStep) + partSys.s0[pOff + PART_XPOS];
    partSys.sM[pOff + PART_YPOS] = (partSys.sdot[pOff + PART_YPOS] / timeStep) + partSys.s0[pOff + PART_YPOS];
    partSys.sM[pOff + PART_ZPOS] = (partSys.sdot[pOff + PART_ZPOS] / timeStep) + partSys.s0[pOff + PART_ZPOS];
    partSys.sM[pOff + PART_XVEL] = (partSys.sdot[pOff + PART_XVEL] / timeStep) + partSys.s0[pOff + PART_XVEL];
    partSys.sM[pOff + PART_YVEL] = (partSys.sdot[pOff + PART_YVEL] / timeStep) + partSys.s0[pOff + PART_YVEL];
    partSys.sM[pOff + PART_ZVEL] = (partSys.sdot[pOff + PART_ZVEL] / timeStep) + partSys.s0[pOff + PART_ZVEL];
    partSys.sM[pOff + PART_R] = (partSys.sdot[pOff + PART_R] / timeStep) + partSys.s0[pOff + PART_R];
    partSys.sM[pOff + PART_G] = (partSys.sdot[pOff + PART_G] / timeStep) + partSys.s0[pOff + PART_G];
    partSys.sM[pOff + PART_B] = (partSys.sdot[pOff + PART_B] / timeStep) + partSys.s0[pOff + PART_B];
    partSys.sM[pOff + PART_R_VEL] = (partSys.sdot[pOff + PART_R_VEL] / timeStep) + partSys.s0[pOff + PART_R_VEL];
    partSys.sM[pOff + PART_G_VEL] = (partSys.sdot[pOff + PART_G_VEL] / timeStep) + partSys.s0[pOff + PART_G_VEL];
    partSys.sM[pOff + PART_B_VEL] = (partSys.sdot[pOff + PART_B_VEL] / timeStep) + partSys.s0[pOff + PART_B_VEL];
    partSys.sM[pOff + PART_MASS] = (partSys.sdot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS];
    partSys.sM[pOff + PART_AGE] = partSys.sdot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];

}

VBObox7.prototype.applyMidForces = function (mygl, i, partSys) {

    var rad = 0;

    var pOff = i * PART_MAXVAR;
    partSys.sM[pOff + PART_X_FTOT] = 0;
    partSys.sM[pOff + PART_Y_FTOT] = 0;
    partSys.sM[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.sM);
    partSys.cforcer[F_GRAV_E](pOff, partSys.sM);
    partSys.cforcer[F_GRAV_P](pOff, partSys.sM);
    partSys.cforcer[F_WIND](pOff, partSys.sM);
    partSys.cforcer[F_BUBBLE](pOff, partSys.sM);
    partSys.cforcer[F_DRAG](pOff, partSys.sM);

    partSys.cforcer[F_SPRING](pOff, partSys.sM);
    partSys.cforcer[F_SPH](pOff, partSys.sM);
    partSys.cforcer[F_CHARGE](pOff, partSys.sM);
    partSys.cforcer[F_FIRE](pOff, partSys.sM);

    /*const F_NONE      0       // Non-existent force: ignore this CForcer object
const F_MOUSE     1       // Spring-like connection to the mouse cursor; lets
                // you 'grab' and 'wiggle' one particle(or several).
const F_GRAV_E    2       // Earth-gravity: pulls all particles 'downward'.
const F_GRAV_P    3       // Planetary-gravity; particle-pair (e0,e1) attract
                // each other with force== grav* mass0*mass1/ dist^2
const F_WIND      4       // Blowing-wind-like force-field;fcn of 3D position
const F_BUBBLE    5       // Constant inward force towards centerpoint if
                // particle is > max_radius away from centerpoint.
const F_DRAG      6       // Viscous drag -- proportional to neg. velocity.
const F_SPRING    7       // ties together 2 particles; distance sets force
const F_SPH 8       // a big collection of identical springs; lets you
                // make cloth & rubbery shapes as one force-making
                // object, instead of many many F_SPRING objects.
const F_CHARGE    9       // attract/repel by charge and inverse distance;
                // applies to all charged particles.
const F_MAXKINDS  10      // 'max' is always the LAST name in our list;
*/

}

VBObox7.prototype.midDotFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s2dot[pOff + PART_XPOS] = partSys.sM[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
    partSys.s2dot[pOff + PART_YPOS] = partSys.sM[pOff + PART_YVEL];
    partSys.s2dot[pOff + PART_ZPOS] = partSys.sM[pOff + PART_ZVEL];
    partSys.s2dot[pOff + PART_XVEL] = partSys.sM[pOff + PART_X_FTOT] * partSys.sM[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_YVEL] = partSys.sM[pOff + PART_Y_FTOT] * partSys.sM[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_ZVEL] = partSys.sM[pOff + PART_Z_FTOT] * partSys.sM[pOff + PART_MASS];
    partSys.s2dot[pOff + PART_X_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_Y_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_Z_FTOT] = 0.0;
    partSys.s2dot[pOff + PART_R] = 0;
    partSys.s2dot[pOff + PART_G] = partSys.sM[pOff + PART_G_VEL];
    partSys.s2dot[pOff + PART_B] = partSys.sM[pOff + PART_B_VEL];
    partSys.s2dot[pOff + PART_MASS] = partSys.sM[pOff + PART_MASS_VEL];
    partSys.s2dot[pOff + PART_DIAM] = 0;
    partSys.s2dot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.s2dot[pOff + PART_AGE] = .08;
    //		partSys.s2dot[pOff + PART_AGE]      = (i/partSys.partCount) * 10;  // # of frame-times since creation/initialization
    partSys.s2dot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.s2dot[pOff + PART_MASS_VEL] = partSys.sM[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
    partSys.s2dot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.s2dot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.s2dot[pOff + PART_G_VEL] = partSys.sM[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
    partSys.s2dot[pOff + PART_B_VEL] = partSys.sM[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
    partSys.s2dot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.s2dot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.s2dot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu

}


VBObox7.prototype.render = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s1[pOff + PART_XPOS] = (partSys.s2dot[pOff + PART_XPOS] / timeStep) + partSys.s0[pOff + PART_XPOS];
    partSys.s1[pOff + PART_YPOS] = (partSys.s2dot[pOff + PART_YPOS] / timeStep) + partSys.s0[pOff + PART_YPOS];
    partSys.s1[pOff + PART_ZPOS] = (partSys.s2dot[pOff + PART_ZPOS] / timeStep) + partSys.s0[pOff + PART_ZPOS];
    partSys.s1[pOff + PART_XVEL] = (partSys.s2dot[pOff + PART_XVEL] / timeStep) + partSys.s0[pOff + PART_XVEL];
    partSys.s1[pOff + PART_YVEL] = (partSys.s2dot[pOff + PART_YVEL] / timeStep) + partSys.s0[pOff + PART_YVEL];
    partSys.s1[pOff + PART_ZVEL] = (partSys.s2dot[pOff + PART_ZVEL] / timeStep) + partSys.s0[pOff + PART_ZVEL];
    partSys.s1[pOff + PART_R] = (partSys.s2dot[pOff + PART_R] / timeStep) + partSys.s0[pOff + PART_R];
    partSys.s1[pOff + PART_G] = (partSys.s2dot[pOff + PART_G] / timeStep) + partSys.s0[pOff + PART_G];
    partSys.s1[pOff + PART_B] = (partSys.s2dot[pOff + PART_B] / timeStep) + partSys.s0[pOff + PART_B];
    partSys.s1[pOff + PART_R_VEL] = (partSys.s2dot[pOff + PART_R_VEL] / timeStep) + partSys.s0[pOff + PART_R_VEL];
    partSys.s1[pOff + PART_G_VEL] = (partSys.s2dot[pOff + PART_G_VEL] / timeStep) + partSys.s0[pOff + PART_G_VEL];
    partSys.s1[pOff + PART_B_VEL] = (partSys.s2dot[pOff + PART_B_VEL] / timeStep) + partSys.s0[pOff + PART_B_VEL];
    partSys.s1[pOff + PART_MASS] = (partSys.s2dot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS];
    partSys.s1[pOff + PART_AGE] = partSys.s2dot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];

    //			if(pOff == 0) console.log(partSys.s1[pOff + PART_XPOS]);
    //		console.log(partSys.s1[pOff + PART_ZVEL]);
}

VBObox7.prototype.applyForces3 = function (mygl, i, partSys) {

    var rad = 0;
    var pOff = i * PART_MAXVAR;
    partSys.s1[pOff + PART_X_FTOT] = 0;
    partSys.s1[pOff + PART_Y_FTOT] = 0;
    partSys.s1[pOff + PART_Z_FTOT] = 0;
    partSys.cforcer[F_MOUSE](pOff, partSys.s1);
    partSys.cforcer[F_GRAV_E](pOff, partSys.s1);
    partSys.cforcer[F_GRAV_P](pOff, partSys.s1);
    partSys.cforcer[F_WIND](pOff, partSys.s1);
    partSys.cforcer[F_BUBBLE](pOff, partSys.s1);
    partSys.cforcer[F_DRAG](pOff, partSys.s1);
    partSys.cforcer[F_SPRING](pOff, partSys.s1);
    partSys.cforcer[F_SPH](pOff, partSys.s1);
    partSys.cforcer[F_CHARGE](pOff, partSys.s1);
    partSys.cforcer[F_FIRE](pOff, partSys.s1);

    /*const F_NONE      0       // Non-existent force: ignore this CForcer object
const F_MOUSE     1       // Spring-like connection to the mouse cursor; lets
                // you 'grab' and 'wiggle' one particle(or several).
const F_GRAV_E    2       // Earth-gravity: pulls all particles 'downward'.
const F_GRAV_P    3       // Planetary-gravity; particle-pair (e0,e1) attract
                // each other with force== grav* mass0*mass1/ dist^2
const F_WIND      4       // Blowing-wind-like force-field;fcn of 3D position
const F_BUBBLE    5       // Constant inward force towards centerpoint if
                // particle is > max_radius away from centerpoint.
const F_DRAG      6       // Viscous drag -- proportional to neg. velocity.
const F_SPRING    7       // ties together 2 particles; distance sets force
const F_SPH 8       // a big collection of identical springs; lets you
                // make cloth & rubbery shapes as one force-making
                // object, instead of many many F_SPRING objects.
const F_CHARGE    9       // attract/repel by charge and inverse distance;
                // applies to all charged particles.
const F_MAXKINDS  10      // 'max' is always the LAST name in our list;
*/

}

VBObox7.prototype.revDotFinder = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.sdot[pOff + PART_XPOS] = -partSys.s1[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
    partSys.sdot[pOff + PART_YPOS] = -partSys.s1[pOff + PART_YVEL];
    partSys.sdot[pOff + PART_ZPOS] = -partSys.s1[pOff + PART_ZVEL];
    //			console.log(partSys.s0[pOff + PART_X_FTOT]);
    partSys.sdot[pOff + PART_XVEL] = -partSys.s0[pOff + PART_X_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.sdot[pOff + PART_YVEL] = -partSys.s0[pOff + PART_Y_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.sdot[pOff + PART_ZVEL] = -partSys.s0[pOff + PART_Z_FTOT] * partSys.s1[pOff + PART_MASS];
    partSys.sdot[pOff + PART_X_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Y_FTOT] = 0.0;
    partSys.sdot[pOff + PART_Z_FTOT] = 0.0;
    partSys.sdot[pOff + PART_R] = 0;
    partSys.sdot[pOff + PART_G] = -partSys.s0[pOff + PART_G_VEL];
    partSys.sdot[pOff + PART_B] = -partSys.s0[pOff + PART_B_VEL];
    partSys.sdot[pOff + PART_MASS] = -partSys.s0[pOff + PART_MASS_VEL];
    partSys.sdot[pOff + PART_DIAM] = 0;
    partSys.sdot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
    partSys.sdot[pOff + PART_AGE] = -.08;
    partSys.sdot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
    partSys.sdot[pOff + PART_MASS_VEL] = -partSys.s0[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
    partSys.sdot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
    partSys.sdot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
    partSys.sdot[pOff + PART_G_VEL] = -partSys.s0[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
    partSys.sdot[pOff + PART_B_VEL] = -partSys.s0[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
    partSys.sdot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
    partSys.sdot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
    partSys.sdot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
}

VBObox7.prototype.render2 = function (mygl, i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.s1[pOff + PART_XPOS] = (partSys.s2dot[pOff + PART_XPOS] / timeStep) + partSys.s0[pOff + PART_XPOS] - (((partSys.s1[pOff + PART_XPOS] + (partSys.sdot[pOff + PART_XPOS] / timeStep)) - partSys.s0[pOff + PART_XPOS]) / 2);
    partSys.s1[pOff + PART_YPOS] = (partSys.s2dot[pOff + PART_YPOS] / timeStep) + partSys.s0[pOff + PART_YPOS] - (((partSys.s1[pOff + PART_YPOS] + (partSys.sdot[pOff + PART_YPOS] / timeStep)) - partSys.s0[pOff + PART_YPOS]) / 2);
    partSys.s1[pOff + PART_ZPOS] = (partSys.s2dot[pOff + PART_ZPOS] / timeStep) + partSys.s0[pOff + PART_ZPOS] - (((partSys.s1[pOff + PART_ZPOS] + (partSys.sdot[pOff + PART_ZPOS] / timeStep)) - partSys.s0[pOff + PART_ZPOS]) / 2);
    partSys.s1[pOff + PART_XVEL] = (partSys.s2dot[pOff + PART_XVEL] / timeStep) + partSys.s0[pOff + PART_XVEL] - (((partSys.s1[pOff + PART_XVEL] + (partSys.sdot[pOff + PART_XVEL] / timeStep)) - partSys.s0[pOff + PART_XVEL]) / 2);
    partSys.s1[pOff + PART_YVEL] = (partSys.s2dot[pOff + PART_YVEL] / timeStep) + partSys.s0[pOff + PART_YVEL] - (((partSys.s1[pOff + PART_YVEL] + (partSys.sdot[pOff + PART_YVEL] / timeStep)) - partSys.s0[pOff + PART_YVEL]) / 2);
    partSys.s1[pOff + PART_ZVEL] = (partSys.s2dot[pOff + PART_ZVEL] / timeStep) + partSys.s0[pOff + PART_ZVEL] - (((partSys.s1[pOff + PART_ZVEL] + (partSys.sdot[pOff + PART_ZVEL] / timeStep)) - partSys.s0[pOff + PART_ZVEL]) / 2);
    partSys.s1[pOff + PART_R] = (partSys.s2dot[pOff + PART_R] / timeStep) + partSys.s0[pOff + PART_R] - (((partSys.s1[pOff + PART_R] + (partSys.sdot[pOff + PART_R] / timeStep)) - partSys.s0[pOff + PART_R]) / 2);
    partSys.s1[pOff + PART_G] = (partSys.s2dot[pOff + PART_G] / timeStep) + partSys.s0[pOff + PART_G] - (((partSys.s1[pOff + PART_G] + (partSys.sdot[pOff + PART_G] / timeStep)) - partSys.s0[pOff + PART_G]) / 2);
    partSys.s1[pOff + PART_B] = (partSys.s2dot[pOff + PART_B] / timeStep) + partSys.s0[pOff + PART_B] - (((partSys.s1[pOff + PART_B] + (partSys.sdot[pOff + PART_B] / timeStep)) - partSys.s0[pOff + PART_B]) / 2);
    partSys.s1[pOff + PART_R_VEL] = (partSys.s2dot[pOff + PART_R_VEL] / timeStep) + partSys.s0[pOff + PART_R_VEL] - (((partSys.s1[pOff + PART_R_VEL] + (partSys.sdot[pOff + PART_R_VEL] / timeStep)) - partSys.s0[pOff + PART_R_VEL]) / 2);
    partSys.s1[pOff + PART_G_VEL] = (partSys.s2dot[pOff + PART_G_VEL] / timeStep) + partSys.s0[pOff + PART_G_VEL] - (((partSys.s1[pOff + PART_G_VEL] + (partSys.sdot[pOff + PART_G_VEL] / timeStep)) - partSys.s0[pOff + PART_G_VEL]) / 2);
    partSys.s1[pOff + PART_B_VEL] = (partSys.s2dot[pOff + PART_B_VEL] / timeStep) + partSys.s0[pOff + PART_B_VEL] - (((partSys.s1[pOff + PART_B_VEL] + (partSys.sdot[pOff + PART_B_VEL] / timeStep)) - partSys.s0[pOff + PART_B_VEL]) / 2);
    partSys.s1[pOff + PART_MASS] = (partSys.s2dot[pOff + PART_MASS] / timeStep) + partSys.s0[pOff + PART_MASS] - (((partSys.s1[pOff + PART_MASS] + (partSys.sdot[pOff + PART_MASS] / timeStep)) - partSys.s0[pOff + PART_MASS]) / 2);
    partSys.s1[pOff + PART_AGE] = partSys.s2dot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];

    //			if(pOff == 0) console.log(partSys.s1[pOff + PART_XPOS]);
    //		console.log(partSys.s1[pOff + PART_ZVEL]);
}

VBObox7.prototype.PartSys_constrain1 = function (i, partSys) {
    var pOff = i * PART_MAXVAR;
    partSys.c0[0](i, partSys.s1, partSys.s0);
    partSys.c0[1](pOff, partSys.s1);
}

VBObox7.prototype.PartSys_render = function (myGL, partSys) {

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
    // 	myGL.uniform1i(this.u_runModeID, myRunMode);		// run/step/pause the particle system

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
    myGL.bufferSubData(myGL.ARRAY_BUFFER, 0, partSys.s0);
    myGL.drawArrays(myGL.POINTS, 0, partSys.partCount);

}


VBObox7.prototype.reset = function (partSys) {
    /*for (i = 0; i < partSys.partCount; i++){
                if(s0[i* PART_MAXVAR + PART_XVEL] > 0.0) s0[i* PART_MAXVAR + PART_XVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_XVEL] -= INIT_VEL;
                if(s0[i* PART_MAXVAR + PART_YVEL] > 0.0) s0[i* PART_MAXVAR + PART_YVEL] += INIT_VEL*.8; else s0[i* PART_MAXVAR + PART_YVEL] -= INIT_VEL*.8;
                if(s0[i* PART_MAXVAR + PART_ZVEL] > 0.0) s0[i* PART_MAXVAR + PART_ZVEL] += INIT_VEL; else s0[i* PART_MAXVAR + PART_ZVEL] -= INIT_VEL;
                }*/
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s0[pOff + PART_XVEL] > 0) {
            partSys.s0[pOff + PART_XVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_XVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_YVEL] > 0) {
            partSys.s0[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_YVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_ZVEL] > 0) {
            partSys.s0[pOff + PART_ZVEL] += (0.2 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_ZVEL] -= (0.2 + 0.8 * Math.random()) * INIT_VEL;
    }
}

VBObox7.prototype.resetWind = function (partSys) {
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        if (partSys.s0[pOff + PART_XVEL] > 0) {
            partSys.s0[pOff + PART_XVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_XVEL] -= (0.3 + 0.8 * Math.random()) * INIT_VEL;
        if (partSys.s0[pOff + PART_YVEL] > 0) {
            partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
        } else partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8 * Math.random()) * INIT_VEL;
    }
}

VBObox7.prototype.resetSpring = function (partSys) {
    for (var i = 1; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        partSys.s0[pOff + PART_ZPOS] -= springRun * (.002 * Math.pow(i, 2));
        /*			}
                else partSys.s0[pOff + PART_XVEL] -= (0.3 + 0.8*Math.random())*INIT_VEL;
                if(  partSys.s0[pOff + PART_YVEL] > 0) {
                       partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8*Math.random())*INIT_VEL;
                    }
                else partSys.s0[pOff + PART_YVEL] += (0.3 + 0.8*Math.random())*INIT_VEL;
            }*/
    }
}