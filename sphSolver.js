function FluidVBObox() {
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

    this.numVox = 2.0/h;


};


FluidVBObox.prototype.init = function (myGL) {
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

FluidVBObox.prototype.adjust1 = function (mygl, partSys) {
    mygl.useProgram(this.shaderLoc);	// In the GPU, SELECT our already-compiled
    partSys.minDensity = 1000000;
    partSys.maxDensity = 1;
    this.applyForces(mygl, partSys);
    this.s0dotFinder(mygl, partSys);
    this.findFirstS1(mygl, partSys);
    this.applyS1Forces(mygl, partSys);
    this.s1DotFinder(mygl, partSys);
    this.findBackS2(mygl, partSys);
    this.render(mygl, partSys);
    this.PartSys_constrain1(partSys);
    partSys.buffer = partSys.s0.slice();
    partSys.s0 = partSys.s1.slice();
    partSys.s1 = partSys.buffer.slice();
    this.vertices = partSys.s0;
    this.PartSys_render(gl, partSys);

}
FluidVBObox.prototype.applyForces = function (mygl, partSys) {

        partSys.cforcer[F_MOUSE](partSys.partCount, partSys.s0);
        partSys.cforcer[F_GRAV_E](partSys.partCount, partSys.s0);
        partSys.cforcer[F_GRAV_P](partSys.partCount, partSys.s0);
        partSys.cforcer[F_WIND](partSys.partCount, partSys.s0);
        partSys.cforcer[F_BUBBLE](partSys.partCount, partSys.s0);
        partSys.cforcer[F_DRAG](partSys.partCount, partSys.s0);
        partSys.cforcer[F_SPRING](partSys.partCount, partSys.s0);
        let results = partSys.cforcer[F_SPH](partSys.partCount, partSys.s0);
        partSys.cforcer[F_CHARGE](partSys.partCount, partSys.s0);
        partSys.cforcer[F_FIRE](partSys.partCount, partSys.s0);

        partSys.minDensity = results[0];
        partSys.maxDensity = results[1];

}


FluidVBObox.prototype.s0dotFinder = function (mygl, partSys) {

    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;
  //      console.log("s0s: ");
    //    console.log("x: " + partSys.s0[pOff + PART_XPOS] + "\ny: " + partSys.s0[pOff + PART_YPOS] + "\nz: " +
      //      partSys.s0[pOff + PART_ZPOS]);
        partSys.s0dot[pOff + PART_XPOS] = partSys.s0[pOff + PART_XVEL];		// 0.0 <= randomRound() < 1.0
        partSys.s0dot[pOff + PART_YPOS] = partSys.s0[pOff + PART_YVEL];
        partSys.s0dot[pOff + PART_ZPOS] = partSys.s0[pOff + PART_ZVEL];
        partSys.s0dot[pOff + PART_XVEL] = partSys.s0[pOff + PART_X_FTOT] * partSys.s0[pOff + PART_MASS];
        partSys.s0dot[pOff + PART_YVEL] = partSys.s0[pOff + PART_Y_FTOT] * partSys.s0[pOff + PART_MASS];
        partSys.s0dot[pOff + PART_ZVEL] = partSys.s0[pOff + PART_Z_FTOT] * partSys.s0[pOff + PART_MASS];
        partSys.s0dot[pOff + PART_X_FTOT] = 0.0;
        partSys.s0dot[pOff + PART_Y_FTOT] = 0.0;
        partSys.s0dot[pOff + PART_Z_FTOT] = 0.0;
        partSys.s0dot[pOff + PART_R] = 0;
        partSys.s0dot[pOff + PART_G] = partSys.s0[pOff + PART_G_VEL];
        partSys.s0dot[pOff + PART_B] = partSys.s0[pOff + PART_B_VEL];
        partSys.s0dot[pOff + PART_MASS] = partSys.s0[pOff + PART_MASS_VEL];
        partSys.s0dot[pOff + PART_DIAM] = 0;
        partSys.s0dot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
        partSys.s0dot[pOff + PART_AGE] = .08;
        partSys.s0dot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
        partSys.s0dot[pOff + PART_MASS_VEL] = partSys.s0[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
        partSys.s0dot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
        partSys.s0dot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
        partSys.s0dot[pOff + PART_G_VEL] = partSys.s0[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
        partSys.s0dot[pOff + PART_B_VEL] = partSys.s0[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
        partSys.s0dot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
        partSys.s0dot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
        partSys.s0dot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu

 //       console.log("Dots: ");
   //     console.log("x: " + partSys.s0dot[pOff + PART_XPOS] + "\ny: " + partSys.s0dot[pOff + PART_YPOS] + "\nz: " +
     //       partSys.s0dot[pOff + PART_ZPOS]);
    }

}

FluidVBObox.prototype.findFirstS1 = function (mygl, partSys) {

    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;
        partSys.s1[pOff + PART_XPOS] = (partSys.s0dot[pOff + PART_XPOS] * timeStep) + partSys.s0[pOff + PART_XPOS];
        partSys.s1[pOff + PART_YPOS] = (partSys.s0dot[pOff + PART_YPOS] * timeStep) + partSys.s0[pOff + PART_YPOS];
        partSys.s1[pOff + PART_ZPOS] = (partSys.s0dot[pOff + PART_ZPOS] * timeStep) + partSys.s0[pOff + PART_ZPOS];
        partSys.s1[pOff + PART_XVEL] = (partSys.s0dot[pOff + PART_XVEL] * timeStep) + partSys.s0[pOff + PART_XVEL];
        partSys.s1[pOff + PART_YVEL] = (partSys.s0dot[pOff + PART_YVEL] * timeStep) + partSys.s0[pOff + PART_YVEL];
        partSys.s1[pOff + PART_ZVEL] = (partSys.s0dot[pOff + PART_ZVEL] * timeStep) + partSys.s0[pOff + PART_ZVEL];
        partSys.s1[pOff + PART_R] = (partSys.s0dot[pOff + PART_R] * timeStep) + partSys.s0[pOff + PART_R];
        partSys.s1[pOff + PART_G] = (partSys.s0dot[pOff + PART_G] * timeStep) + partSys.s0[pOff + PART_G];
        partSys.s1[pOff + PART_B] = (partSys.s0dot[pOff + PART_B] * timeStep) + partSys.s0[pOff + PART_B];
        partSys.s1[pOff + PART_R_VEL] = (partSys.s0dot[pOff + PART_R_VEL] * timeStep) + partSys.s0[pOff + PART_R_VEL];
        partSys.s1[pOff + PART_G_VEL] = (partSys.s0dot[pOff + PART_G_VEL] * timeStep) + partSys.s0[pOff + PART_G_VEL];
        partSys.s1[pOff + PART_B_VEL] = (partSys.s0dot[pOff + PART_B_VEL] * timeStep) + partSys.s0[pOff + PART_B_VEL];
        partSys.s1[pOff + PART_MASS] = (partSys.s0dot[pOff + PART_MASS] * timeStep) + partSys.s0[pOff + PART_MASS];
        partSys.s1[pOff + PART_AGE] = partSys.s0dot[pOff + PART_AGE] + partSys.s0[pOff + PART_AGE];
 //       console.log("first s1: ");
   //     console.log("x: " + partSys.s1[pOff + PART_XPOS] + "\ny: " + partSys.s1[pOff + PART_YPOS] + "\nz: " +
     //       partSys.s1[pOff + PART_ZPOS]);

    }
}

FluidVBObox.prototype.applyS1Forces = function (mygl, partSys) {

        partSys.cforcer[F_MOUSE](partSys.partCount, partSys.s1);
        partSys.cforcer[F_GRAV_E](partSys.partCount, partSys.s1);
        partSys.cforcer[F_GRAV_P](partSys.partCount, partSys.s1);
        partSys.cforcer[F_WIND](partSys.partCount, partSys.s1);
        partSys.cforcer[F_BUBBLE](partSys.partCount, partSys.s1);
        partSys.cforcer[F_DRAG](partSys.partCount, partSys.s1);
        partSys.cforcer[F_SPRING](partSys.partCount, partSys.s1);
        let results = partSys.cforcer[F_SPH](partSys.partCount, partSys.s1);
        partSys.cforcer[F_CHARGE](partSys.partCount, partSys.s1);
        partSys.cforcer[F_FIRE](partSys.partCount, partSys.s1);

    partSys.minDensity = results[0];
    partSys.maxDensity = results[1];
}

FluidVBObox.prototype.s1DotFinder = function (mygl, partSys) {

    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;
        partSys.s1dot[pOff + PART_XPOS] = partSys.s1[pOff + PART_XVEL];
        partSys.s1dot[pOff + PART_YPOS] = partSys.s1[pOff + PART_YVEL];
        partSys.s1dot[pOff + PART_ZPOS] = partSys.s1[pOff + PART_ZVEL];
        partSys.s1dot[pOff + PART_XVEL] = partSys.s1[pOff + PART_X_FTOT] * partSys.s1[pOff + PART_MASS];
        partSys.s1dot[pOff + PART_YVEL] = partSys.s1[pOff + PART_Y_FTOT] * partSys.s1[pOff + PART_MASS];
        partSys.s1dot[pOff + PART_ZVEL] = partSys.s1[pOff + PART_Z_FTOT] * partSys.s1[pOff + PART_MASS];
        partSys.s1dot[pOff + PART_X_FTOT] = 0.0;
        partSys.s1dot[pOff + PART_Y_FTOT] = 0.0;
        partSys.s1dot[pOff + PART_Z_FTOT] = 0.0;
        partSys.s1dot[pOff + PART_R] = 0;
        partSys.s1dot[pOff + PART_G] = partSys.s1[pOff + PART_G_VEL];
        partSys.s1dot[pOff + PART_B] = partSys.s1[pOff + PART_B_VEL];
        partSys.s1dot[pOff + PART_MASS] = partSys.s1[pOff + PART_MASS_VEL];
        partSys.s1dot[pOff + PART_DIAM] = 0;
        partSys.s1dot[pOff + PART_RENDMODE] = 0; // 0,1,2 or 3
        partSys.s1dot[pOff + PART_AGE] = .08;
        partSys.s1dot[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
        partSys.s1dot[pOff + PART_MASS_VEL] = partSys.s1[pOff + PART_MASS_FTOT];  // time-rate-of-change of mass.
        partSys.s1dot[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
        partSys.s1dot[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
        partSys.s1dot[pOff + PART_G_VEL] = partSys.s1[pOff + PART_G_FTOT];  // time-rate-of-change of color:grn
        partSys.s1dot[pOff + PART_B_VEL] = partSys.s1[pOff + PART_B_FTOT]; // time-rate-of-change of color:blu
        partSys.s1dot[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
        partSys.s1dot[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
        partSys.s1dot[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
 //       console.log("1Dots: ");
   //     console.log("x: " + partSys.s1dot[pOff + PART_XPOS] + "\ny: " + partSys.s1dot[pOff + PART_YPOS] + "\nz: " +
     //       partSys.s1dot[pOff + PART_ZPOS]);
    }

}

FluidVBObox.prototype.findBackS2 = function (mygl, partSys) {
    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;
        partSys.s2[pOff + PART_XPOS] = -(partSys.s1dot[pOff + PART_XPOS] * timeStep) + partSys.s1[pOff + PART_XPOS];
        partSys.s2[pOff + PART_YPOS] = -(partSys.s1dot[pOff + PART_YPOS] * timeStep) + partSys.s1[pOff + PART_YPOS];
        partSys.s2[pOff + PART_ZPOS] = -(partSys.s1dot[pOff + PART_ZPOS] * timeStep) + partSys.s1[pOff + PART_ZPOS];
        partSys.s2[pOff + PART_XVEL] = -(partSys.s1dot[pOff + PART_XVEL] * timeStep) + partSys.s1[pOff + PART_XVEL];
        partSys.s2[pOff + PART_YVEL] = -(partSys.s1dot[pOff + PART_YVEL] * timeStep) + partSys.s1[pOff + PART_YVEL];
        partSys.s2[pOff + PART_ZVEL] = -(partSys.s1dot[pOff + PART_ZVEL] * timeStep) + partSys.s1[pOff + PART_ZVEL];
        partSys.s2[pOff + PART_R] = -(partSys.s1dot[pOff + PART_R] * timeStep) + partSys.s1[pOff + PART_R];
        partSys.s2[pOff + PART_G] = -(partSys.s1dot[pOff + PART_G] * timeStep) + partSys.s1[pOff + PART_G];
        partSys.s2[pOff + PART_B] = -(partSys.s1dot[pOff + PART_B] * timeStep) + partSys.s1[pOff + PART_B];
        partSys.s2[pOff + PART_R_VEL] = -(partSys.s1dot[pOff + PART_R_VEL] * timeStep) + partSys.s1[pOff + PART_R_VEL];
        partSys.s2[pOff + PART_G_VEL] = -(partSys.s1dot[pOff + PART_G_VEL] * timeStep) + partSys.s1[pOff + PART_G_VEL];
        partSys.s2[pOff + PART_B_VEL] = -(partSys.s1dot[pOff + PART_B_VEL] * timeStep) + partSys.s1[pOff + PART_B_VEL];
        partSys.s2[pOff + PART_MASS] = -(partSys.s1dot[pOff + PART_MASS] * timeStep) + partSys.s1[pOff + PART_MASS];
        partSys.s2[pOff + PART_AGE] = -partSys.s1dot[pOff + PART_AGE] + partSys.s1[pOff + PART_AGE];
     //   console.log("s2: ");
       // console.log("x: " + partSys.s2[pOff + PART_XPOS] + "\ny: " + partSys.s2[pOff + PART_YPOS] + "\nz: " +
         //   partSys.s2[pOff + PART_ZPOS]);
    }
}

FluidVBObox.prototype.render = function (mygl, partSys) {

    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;
        partSys.s1[pOff + PART_XPOS] += (partSys.s0[pOff + PART_XPOS] - partSys.s2[pOff + PART_XPOS]) / 2;
        partSys.s1[pOff + PART_YPOS] += (partSys.s0[pOff + PART_YPOS] - partSys.s2[pOff + PART_YPOS]) / 2;
        partSys.s1[pOff + PART_ZPOS] += (partSys.s0[pOff + PART_ZPOS] - partSys.s2[pOff + PART_ZPOS]) / 2;
        partSys.s1[pOff + PART_XVEL] += (partSys.s0[pOff + PART_XVEL] - partSys.s2[pOff + PART_XVEL]) / 2;
        partSys.s1[pOff + PART_YVEL] += (partSys.s0[pOff + PART_YVEL] - partSys.s2[pOff + PART_YVEL]) / 2;
        partSys.s1[pOff + PART_ZVEL] += (partSys.s0[pOff + PART_ZVEL] - partSys.s2[pOff + PART_ZVEL]) / 2;
        partSys.s1[pOff + PART_R] += (partSys.s0[pOff + PART_R] - partSys.s2[pOff + PART_R]) / 2;
        partSys.s1[pOff + PART_G] += (partSys.s0[pOff + PART_G] - partSys.s2[pOff + PART_G]) / 2;
        partSys.s1[pOff + PART_B] += (partSys.s0[pOff + PART_B] - partSys.s2[pOff + PART_B]) / 2;
        partSys.s1[pOff + PART_R_VEL] += (partSys.s0[pOff + PART_R_VEL] - partSys.s2[pOff + PART_R_VEL]) / 2;
        partSys.s1[pOff + PART_G_VEL] += (partSys.s0[pOff + PART_G_VEL] - partSys.s2[pOff + PART_G_VEL]) / 2;
        partSys.s1[pOff + PART_B_VEL] += (partSys.s0[pOff + PART_B_VEL] - partSys.s2[pOff + PART_B_VEL]) / 2;
        partSys.s1[pOff + PART_MASS] += (partSys.s0[pOff + PART_MASS] - partSys.s2[pOff + PART_MASS]) / 2;
        partSys.s1[pOff + PART_AGE] += (partSys.s0[pOff + PART_AGE] - partSys.s0dot[pOff + PART_AGE]);
        partSys.s1[pOff + PART_DENSITY] = partSys.s0[pOff + PART_DENSITY];
        partSys.s1[pOff + PART_VOX_X] = Math.max(Math.floor(partSys.s1[pOff + PART_XPOS] / (h * 2)), 0);
        partSys.s1[pOff + PART_VOX_Y] = Math.max(Math.floor(partSys.s1[pOff + PART_YPOS] / (h * 2)), 0);
        partSys.s1[pOff + PART_VOX_Z] = Math.max(Math.floor(partSys.s1[pOff + PART_ZPOS] / (h * 2)), 0);
  //      console.log("1s: ");
    //    console.log("x: " + (partSys.s0[pOff + PART_XPOS] - partSys.s2[pOff + PART_XPOS]) + "\ny: " +
      //      (partSys.s0[pOff + PART_YPOS] - partSys.s2[pOff + PART_YPOS]) + "\nz: " +
        //    (partSys.s0[pOff + PART_ZPOS] - partSys.s2[pOff + PART_ZPOS]));

    }
}

FluidVBObox.prototype.PartSys_constrain1 = function (partSys) {

    for (var i = 0; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;
        partSys.c0[0](i, partSys.s1, partSys.s0);
        partSys.c0[1](pOff, partSys.s1);
    //    console.log("s1's: \n x:" + partSys.s1[pOff + PART_XPOS] + "\ny: " + partSys.s1[pOff + PART_YPOS]
      //      + "\n z:" + partSys.s1[pOff + PART_ZPOS])
    }
}

FluidVBObox.prototype.PartSys_render = function (myGL, partSys) {

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


FluidVBObox.prototype.reset = function (partSys) {

    for (var i = 0; i < partSys.partCount; i++) {
        partSys.s0[pOff + PART_XVEL] += (-0.5 + Math.random()) * .02;
        partSys.s0[pOff + PART_YVEL] += (-0.5 + Math.random()) * .02;
        partSys.s0[pOff + PART_ZVEL] += (-0.5 + Math.random()) * .02;
    }
}
FluidVBObox.prototype.resetWind = function (partSys) {
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
FluidVBObox.prototype.resetSpring = function (partSys) {
    for (var i = 1; i < partSys.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        partSys.s0[pOff + PART_ZPOS] -= springRun * (.002 * Math.pow(i, 2));
    }
}

function createArray(length) {
    var arr = new Array(length || 0),
        i = length;

    if (arguments.length > 1) {
        var args = Array.prototype.slice.call(arguments, 1);
        while(i--) arr[length-1 - i] = createArray.apply(this, args);
    }

    return arr;
}
