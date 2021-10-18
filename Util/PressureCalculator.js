findPressureForce = function (s, partCount) {

    for (let j = 0; j < partCount; j++) {
        let p = j * PART_MAXVAR;
        var firstComponent = [0, 0, 0];
        var secondComponent = [0, 0, 0];
        for (let i = 0; i < partCount; i++) {
            let pOff = i * PART_MAXVAR;
            if (sameVoxel(s, p, pOff)) {
                if (pOff !== p) {
                    let {xDist, yDist, zDist, realDistance} = calculateDistances(s, p, pOff);
                    if (realDistance < 2 * h) {
                        //Kernel function = 45 * (2 - (r)^2)
                        //                    (64 *pi)
                        //let kernelDer = -15 * (Math.pow(h - realDistance, 3))/ (Math.PI * Math.pow(h, 6)); //
                        // let kernelDer = -45 * (2 - Math.pow((realDistance / h), 2)) / (64 * h * Math.PI);
                        let kernelDer = calculateKernel(realDistance)

                        let xdistanceComponent = xDist / realDistance;
                        let ydistanceComponent = yDist / realDistance;
                        let zdistanceComponent = zDist / realDistance;
                        //   if (distance > 0) {

                        let forceOfPressure = -.2 * s[pOff + PART_MASS]/s[pOff + PART_DENSITY] * (s[p + PART_PRESSURE]/Math.pow(s[p + PART_DENSITY], 2)
                            + s[pOff + PART_PRESSURE]/Math.pow(s[pOff + PART_DENSITY], 2));
                        s[p + PART_X_FTOT] -= forceOfPressure * kernelDer * xdistanceComponent;
                        s[p + PART_Y_FTOT] -= forceOfPressure * kernelDer * ydistanceComponent;
                        s[p + PART_Z_FTOT] -= forceOfPressure * kernelDer * zdistanceComponent;
                        /*                        firstComponent[0] += s[pOff + PART_MASS] * kernelDer * (-1 * xdistanceComponent);
                                                secondComponent[0] += s[pOff + PART_MASS] * ((s[pOff + PART_DENSITY] - DENSITY_CONST) / ((s[pOff + PART_DENSITY] * (s[pOff + PART_DENSITY])))) * kernelDer * (-1 * xdistanceComponent);

                                                //   if (distance > 0) {
                                                firstComponent[1] += s[pOff + PART_MASS] * kernelDer * (-1 * ydistanceComponent);
                                                secondComponent[1] += s[pOff + PART_MASS] * ((s[pOff + PART_DENSITY] - DENSITY_CONST) / ((s[pOff + PART_DENSITY] * (s[pOff + PART_DENSITY])))) * kernelDer * (-1 * ydistanceComponent);

                                                //   if (distance > 0) {
                                                firstComponent[2] += s[pOff + PART_MASS] * kernelDer * (-1 * zdistanceComponent);
                                                secondComponent[2] += s[pOff + PART_MASS] * ((s[pOff + PART_DENSITY] - DENSITY_CONST) / ((s[pOff + PART_DENSITY] * (s[pOff + PART_DENSITY])))) * kernelDer * (-1 * zdistanceComponent);
                        */
                        //}
                        /*   else{
                                firstComponent[pos] -= s[pOff + PART_MASS] * kernelDer * distanceComponent;
                              secondComponent[pos] -= s[pOff + PART_MASS] * ((s[pOff + PART_DENSITY] - DENSITY_CONST) / (s[pOff + PART_DENSITY])) * kernelDer * distanceComponent;
                           }*/

                    }
                }
            }
        }
        //       let coefficient = k * s[p + PART_MASS];
        //     let firstDensCo = ((s[p + PART_DENSITY] - DENSITY_CONST) / (s[p + PART_DENSITY]* s[p + PART_DENSITY]));
        //    s[p + PART_X_FTOT] += coefficient * ((firstDensCo * firstComponent[0]) + secondComponent[0]);
        //  s[p + PART_Y_FTOT] += coefficient * ((firstDensCo * firstComponent[1]) + secondComponent[1]);
        //s[p + PART_Z_FTOT] += coefficient * ((firstDensCo * firstComponent[2]) + secondComponent[2]);
        //            console.log('x ' + s[p + PART_X_FTOT] + '\ny ' + s[p + PART_Y_FTOT] + '\nz ' + s[p + PART_Z_FTOT]);
    }


    // Ideas: calculate current pressure and calculate pressure gradient using that
    //   current pressure = speed of sound ^ 2 * (density - goal density)
    // sum all
    //      FPij = - mj/densej * ( pi / densei^2  + pj/densej^2) * kernel

    //current formula  -mass * kernel * (dense - goalDense/ dense^2)
}