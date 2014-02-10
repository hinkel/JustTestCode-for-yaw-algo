#include "board.h"
#include "mw.h"

typedef struct fp_vector
{
    float X, Y, Z;
} t_fp_vector_def;

typedef union
{
    float A[3];
    t_fp_vector_def V;
} t_fp_vector;

float   accSmooth[3], ACC_speed[2], accADC[3], gyroADC[3], magADCfloat[3], angle[2] = { 0, 0 }; // absolute angle inclination in multiple of 0.1 degree    180 deg = 1800
float   BaroAlt, EstAlt, AltHold, vario, ACCDeltaTimeINS = 0;;
int32_t sonarAlt;
int16_t BaroP, BaroI, BaroD;
bool    newbaroalt, GroundAltInitialized, HaveNewMag;

// **************
// gyro+acc IMU
// **************
/*
 * Sensor data rate:
 * Baro  - 25 Hz  - 40 ms  | 50 ms
 * Accel - 400 Hz - 2.5 ms | 10 ms
 * Mag   - 30 Hz  - 4.5 ms | 40 ms
 * Gyro  - 760 Hz - 1.3 ms | 10 ms
 */

float  gyroData[3] = { 0, 0, 0 };                         
static uint8_t SmoothingFactor[3]  = { 0, 0, 0 };
static bool    GyroSmoothing;
static void getEstimatedAttitude(void);

void imuInit(void)                                                            // Initialize & precalculate some values here
{
    if (cfg.gy_smrll || cfg.gy_smptc || cfg.gy_smyw)
    {
        SmoothingFactor[ROLL]  = cfg.gy_smrll;
        SmoothingFactor[PITCH] = cfg.gy_smptc;
        SmoothingFactor[YAW]   = cfg.gy_smyw;
        GyroSmoothing          = true;
    } else GyroSmoothing = false;

#ifdef MAG
    if (sensors(SENSOR_MAG)) Mag_init();
#endif
}

void computeIMU(void)
{
    static  float LastGyroSmooth[3] = { 0.0f, 0.0f, 0.0f };
    static  int16_t triywavg[4];
    static  uint8_t triywavgpIDX = 0;
    uint8_t axis, i;
    float   flttmp;

    if (MpuSpecial) GETMPU6050();
    else
    {
        gyro.temperature(&telemTemperature1);                                 // Read out gyro temperature
        Gyro_getADC();                                                        // Also feeds gyroData
        if (sensors(SENSOR_ACC)) ACC_getADC();
    }
    if(cfg.acc_calibrated) getEstimatedAttitude();                            // acc_calibrated just can turn true if acc present.
    
    if(cfg.mixerConfiguration == MULTITYPE_TRI && cfg.gy_smyw)                // Moving average for yaw in tri mode
    {
        triywavg[triywavgpIDX] = (int16_t)gyroData[YAW]; triywavgpIDX++;
        if (triywavgpIDX == 4) triywavgpIDX = 0;
        flttmp = 0;
        for (i = 0; i < 4; i++) flttmp += triywavg[i];
        gyroData[YAW] = flttmp * 0.25f;
    }

    if (GyroSmoothing)
    {
        for (axis = 0; axis < 3; axis++)
        {
            if (SmoothingFactor[axis] > 1)                                    // Circumvent useless action
            {
                flttmp               = (float)SmoothingFactor[axis];
                gyroData[axis]       = ((LastGyroSmooth[axis] * (flttmp - 1.0f)) + gyroData[axis]) / flttmp;
                LastGyroSmooth[axis] = gyroData[axis];
            }
        }
    }
}

void RotGravAndMag(struct fp_vector *Grav, struct fp_vector *Mag, float *delta)// Rotate vectors according to the gyro delta
{
    struct    fp_vector v_tmp = *Grav;
    float     mat[3][3], cosx, sinx, cosy, siny, cosz, sinz, coszcosx, coszcosy, sinzcosx, coszsinx, sinzsinx;
    cosx      =  cosf(-delta[PITCH]);
    sinx      =  sinf(-delta[PITCH]);
    cosy      =  cosf(delta[ROLL]);
    siny      =  sinf(delta[ROLL]);
    cosz      =  cosf(delta[YAW]);
    sinz      =  sinf(delta[YAW]);
    coszcosx  =  cosz * cosx;
    coszcosy  =  cosz * cosy;
    sinzcosx  =  sinz * cosx;
    coszsinx  =  sinx * cosz;
    sinzsinx  =  sinx * sinz;
    mat[0][0] =  coszcosy;
    mat[0][1] =  sinz * cosy;
    mat[0][2] = -siny;
    mat[1][0] = (coszsinx * siny) - sinzcosx;
    mat[1][1] = (sinzsinx * siny) + (coszcosx);
    mat[1][2] =  cosy * sinx;
    mat[2][0] = (coszcosx * siny) + (sinzsinx);
    mat[2][1] = (sinzcosx * siny) - (coszsinx);
    mat[2][2] =  cosy * cosx;
    Grav->X   =  v_tmp.X * mat[0][0] + v_tmp.Y * mat[1][0] + v_tmp.Z * mat[2][0];
    Grav->Y   =  v_tmp.X * mat[0][1] + v_tmp.Y * mat[1][1] + v_tmp.Z * mat[2][1];
    Grav->Z   =  v_tmp.X * mat[0][2] + v_tmp.Y * mat[1][2] + v_tmp.Z * mat[2][2];
    if (!cfg.mag_calibrated) return;                                                  // mag_calibrated can just be true if MAG present
    v_tmp     = *Mag;                                                                 // Proceed here if mag present and calibrated
    Mag->X    = v_tmp.X * mat[0][0] + v_tmp.Y * mat[1][0] + v_tmp.Z * mat[2][0];      // That saves a butload of recalculating matrix values
    Mag->Y    = v_tmp.X * mat[0][1] + v_tmp.Y * mat[1][1] + v_tmp.Z * mat[2][1];
    Mag->Z    = v_tmp.X * mat[0][2] + v_tmp.Y * mat[1][2] + v_tmp.Z * mat[2][2];
}

// *Somehow* modified by me..
// Besides mwii credit must go to Sebbi, BRM! Hopefully they condone mentioning them above my trash.
// Sebbi for his rotation of the acc vector and BRM for his normalization ideas.
static void getEstimatedAttitude(void)
{
    static t_fp_vector EstG, EstM;
    static float    accLPFALT[3], accLPFGPS[3], Tilt_25deg, INVsens1G;
    static float    INV_GYR_CMPF_FACTOR, INV_GYR_CMPFM_FACTOR, ACC_GPS_RC, ACC_ALT_RC, ACC_RC;
    static uint32_t previousT, UpsDwnTimer;
    static bool     init = false;
    float           tmp[3], scale, deltaGyroAngle[3], ACCALTFac, ACCGPSFac, ACCFac, rollRAD, pitchRAD;
    float           NormFst = 0.0f, NormScnd, NormR, A, B, cr, sr, cp, sp, cy, sy, spcy, spsy, acc_south, acc_west, acc_up, FAC;
    uint8_t         i;
    uint32_t        tmpu32, currentT = micros();
    ACCDeltaTimeINS = (float)(currentT - previousT) * 0.000001f;
    previousT       = currentT;
    if(!init)                                                                 // Setup variables & constants
    {
        init = true;
        INVsens1G            = 1.0f / cfg.sens_1G;
        Tilt_25deg           = cosf(25.0f * RADX);
        INV_GYR_CMPF_FACTOR  = 1.0f / (float)(cfg.gy_cmpf  + 1);              // Default 400
        INV_GYR_CMPFM_FACTOR = 1.0f / (float)(cfg.gy_cmpfm + 1);              // Default 200
        ACC_ALT_RC           = 0.5f / (M_PI * cfg.acc_altlpfhz);              // Default 10 Hz
        ACC_RC               = 0.5f / (M_PI * cfg.acc_lpfhz);                 // Default 0,536 Hz
        ACC_GPS_RC           = 0.5f / (M_PI * cfg.acc_gpslpfhz);              // Default 5 Hz
        for (i = 0; i < 3; i++)                                               // Preset some values to reduce runup time
        {
            tmp[0]       = accADC[i] * INVsens1G;
            accSmooth[i] = tmp[0];
            accLPFGPS[i] = tmp[0];
            accLPFALT[i] = tmp[0];
            EstG.A[i]    = tmp[0] * 0.5f;
        }
    }
    ACCALTFac = ACCDeltaTimeINS / (ACC_ALT_RC + ACCDeltaTimeINS);             // Adjust LPF to cycle time / do Hz cut off
    ACCGPSFac = ACCDeltaTimeINS / (ACC_GPS_RC + ACCDeltaTimeINS);             // Adjust LPF to cycle time / do Hz cut off
    ACCFac    = ACCDeltaTimeINS / (ACC_RC     + ACCDeltaTimeINS);
    scale     = ACCDeltaTimeINS * GyroScale;                                  // SCALE CHANGED TO RAD per SECONDS, DAMN
    for (i = 0; i < 3; i++)
    {
        tmp[0]            = accADC[i]    *  INVsens1G;                        // Reference all to 1G here
        accLPFGPS[i]     += ACCGPSFac    * (tmp[0] - accLPFGPS[i]);           // For GPS
        accLPFALT[i]     += ACCALTFac    * (tmp[0] - accLPFALT[i]);           // For Althold
        accSmooth[i]     += ACCFac       * (tmp[0] - accSmooth[i]);           // For Gyrodrift correction
        NormFst          += accSmooth[i] *  accSmooth[i];
        deltaGyroAngle[i] = gyroADC[i]   *  scale;                            // deltaGyroAngle is in RAD
    }
    RotGravAndMag(&EstG.V, &EstM.V, deltaGyroAngle);                          // Rotate Grav&Mag together to avoid doublecalculation
    NormFst  = sqrtf(NormFst); 
    tmpu32   = (uint32_t)(NormFst * 1000.0f);                                 // Make it "shorter" for comparison in temp variable
    NormScnd = sqrtf(EstG.A[0] * EstG.A[0] + EstG.A[1] * EstG.A[1] + EstG.A[2] * EstG.A[2]);
    if (NormScnd) NormR = 1.0f / NormScnd;
    else          NormR = INVsens1G;                                          // Feed vanilla value in that rare case, to let angle calculation always happen
    for (i = 0; i < 3; i++) tmp[i] = EstG.A[i] * NormR;                       // tmp[0..2] contains normalized EstG
    if (850 < tmpu32 && tmpu32 < 1150)                                        // Gyro drift correction if ACC between 0.85G and 1.15G else skip filter, as EstV already rotated by Gyro
    {
        NormR = 1.0f / NormFst;                                               // We just cmp filter together normalized vectors here. Div 0 not possible here.
        for (i = 0; i < 3; i++) EstG.A[i] = (tmp[i] * (float)cfg.gy_cmpf + accSmooth[i] * NormR) * INV_GYR_CMPF_FACTOR;
    }
    rollRAD      =  atan2f(tmp[0], tmp[2]);                                   // Note: One cycle after successful cmpf is done with old values.
    pitchRAD     =  asinf(-tmp[1]);                                           // Has to have the "wrong sign" relative to angle[PITCH]
    angle[ROLL]  =  rollRAD  * RADtoDEG10;                                    // Use rounded values, eliminate jitter for main PID I and D
    angle[PITCH] = -pitchRAD * RADtoDEG10;
    TiltValue    =  cosf(rollRAD) * cosf(pitchRAD);                           // We do this correctly here
    if (TiltValue >= 0)   UpsDwnTimer = 0;
    else if(!UpsDwnTimer) UpsDwnTimer = currentT + 20000;                     // Use Timer here to make absolutely sure we are upsidedown
    if (UpsDwnTimer && currentT > UpsDwnTimer) UpsideDown = true;
    else UpsideDown = false;
    if (TiltValue > Tilt_25deg) f.SMALL_ANGLES_25 = 1;
    else f.SMALL_ANGLES_25 = 0;
    cr = cosf(rollRAD);
    sr = sinf(rollRAD);
    cp = cosf(pitchRAD);
    sp = sinf(pitchRAD);
    if (cfg.mag_calibrated)                                                   // mag_calibrated can just be true if MAG present
    {
        NormFst = sqrtf(EstM.A[0] * EstM.A[0] + EstM.A[1] * EstM.A[1] + EstM.A[2] * EstM.A[2]);
        if (NormFst) NormR = 1.0f / NormFst;
        else         NormR = 1.0f;                                            // Vanilla value for rare case
        for (i = 0; i < 3; i++) tmp[i] = EstM.A[i] * NormR;                   // tmp[0..2] contains normalized EstM
        if(HaveNewMag)                                                        // Only do Complementary filter when new MAG data are available
        {
            HaveNewMag = false;
            NormFst = sqrtf(magADCfloat[0] * magADCfloat[0] + magADCfloat[1] * magADCfloat[1] + magADCfloat[2] * magADCfloat[2]);
            if (NormFst) NormR = 1.0f / NormFst;
            else         NormR = 1.0f; 
            for (i = 0; i < 3; i++) EstM.A[i] = (tmp[i] * (float)cfg.gy_cmpfm + magADCfloat[i] * NormR) * INV_GYR_CMPFM_FACTOR;
        }
        A  = tmp[1] * cp + tmp[0] * sr * sp + tmp[2] * cr * sp;
        B  = tmp[0] * cr - tmp[2] * sr;
        heading = wrap_180(atan2f(-B, A) * RADtoDEG + magneticDeclination);   // Get rad to Degree and add declination (note: without *10) // Wrap to -180 0 +180 Degree        
    } else heading = 0;                                                       // if no mag or not calibrated do bodyframe below
    tmp[0]    = heading * RADX;                                               // Do GPS INS rotate ACC X/Y to earthframe no centrifugal comp. yet
    cy        = cosf(tmp[0]);
    sy        = sinf(tmp[0]);
    cos_yaw_x = cy;                                                           // Store for general use
    sin_yaw_y = sy;                                                           // Store for general use
    spcy      = sp * cy;
    spsy      = sp * sy;
    FAC       = 980.665f * ACCDeltaTimeINS;                                   // vel factor for normalized output tmp3 = (9.80665f * (float)ACCDeltaTime) / 10000.0f;
    acc_up    = ((-sp) * accLPFALT[1] + sr * cp * accLPFALT[0] + cp * cr * accLPFALT[2]) - 1;// -1G
    if(GroundAltInitialized) vario += acc_up * FAC * constrain(TiltValue + 0.05, 0.5f, 1.0f);// Positive when moving Up. Just do Vario when Baro completely initialized. Empirical hightdrop reduction on tilt.
    acc_south = cp * cy * accLPFGPS[1] + (sr * spcy - cr * sy) * accLPFGPS[0] + ( sr * sy + cr * spcy) * accLPFGPS[2];
    acc_west  = cp * sy * accLPFGPS[1] + (cr * cy + sr * spsy) * accLPFGPS[0] + (-sr * cy + cr * spsy) * accLPFGPS[2];
    ACC_speed[LAT] -= acc_south * FAC;                                        // Positive when moving North cm/sec when no MAG this is speed to the front
    ACC_speed[LON] -= acc_west  * FAC;                                        // Positive when moving East cm/sec when no MAG this is speed to the right
}

#ifdef BARO
///////////////////////////////////////////////
//Crashpilot1000 Mod getEstimatedAltitude ACC//
///////////////////////////////////////////////
#define VarioTabsize 8
void getEstimatedAltitude(void)
{
    static int8_t   VarioTab[VarioTabsize];
    static uint8_t  Vidx = 0, IniStep = 0, IniCnt = 0;
    static uint32_t LastBarotime = 0;
    static float    AvgHz = 0, LastEstAltBaro = 0, SNRcorrect, SNRavg = 0;
    float           NewVal, EstAltBaro;
    uint32_t        TimeTemp;
    uint8_t         i;

    if (!GroundAltInitialized)
    {
        if (newbaroalt)
        {
            TimeTemp     = micros();
            NewVal       = (float)(TimeTemp - LastBarotime);
            LastBarotime = TimeTemp;          
            switch(IniStep)                                                   // Casemachine here for further extension
            {
            case 0:
                IniCnt++;
                if(IniCnt == 50)                                              // Waste 50 Cycles to let things (buffers) settle then ini some vars and proceed
                {
                    for (i = 0; i < VarioTabsize; i++) VarioTab[i] = 0;
                    EstAlt = GroundAlt = vario = 0;
                    IniCnt = SonarStatus = 0;
                    IniStep++;
                }
                break;
            case 1:
                GroundAlt += BaroAlt;
                AvgHz     += NewVal;
                IniCnt++;
                if (IniCnt == 50)                                             // Gather 50 values
                {
                    GroundAlt *= 0.02f;
                    AvgHz      = 50000000.0f / AvgHz;                         // Calculate Average Hz here since we skip Baro temp readout every 2nd read
                    GroundAltInitialized = true;
                }
                break;
            }
        }
    }
    else
    {
        if (SonarStatus) NewVal = sonarAlt;
        switch(SonarStatus)
        {
        case 0:
            SNRavg  = 0;
            IniStep = 0;
            break;
        case 1:
            if (!IniStep)
            {
                IniStep = 1;
                SNRavg  = NewVal;
            }
            else SNRavg += 0.2f * (NewVal - SNRavg);                          // Adjust Average during accepttimer (ca. 550ms so ca. 20 cycles)
            SNRcorrect = EstAlt + GroundAlt - SNRavg;                         // Calculate baro/sonar displacement on 1st contact
            break;
        case 2:
            if (newbaroalt) BaroAlt = (SNRcorrect + NewVal) * cfg.snr_cf + BaroAlt * (1 - cfg.snr_cf); // Set weight / make transition smoother
            break;
        }
        EstAlt += vario * ACCDeltaTimeINS;
        if (newbaroalt)
        {
            EstAltBaro     = BaroAlt - GroundAlt;
            VarioTab[Vidx] = constrain((int16_t)(EstAltBaro - LastEstAltBaro), -127, 127); Vidx++;
            if (Vidx == VarioTabsize) Vidx = 0;
            LastEstAltBaro = EstAltBaro;
            NewVal = 0;
            for (i = 0; i < VarioTabsize; i++) NewVal += (float)VarioTab[i];
            NewVal = (NewVal * AvgHz)/(float)VarioTabsize;
            vario  = vario  * cfg.accz_vcf + NewVal     * (1.0f - cfg.accz_vcf);
            EstAlt = EstAlt * cfg.accz_acf + EstAltBaro * (1.0f - cfg.accz_acf);
            if (cfg.bar_dbg)
            {
                debug[0] = EstAltBaro * 10;
                debug[1] = EstAlt     * 10;
                debug[2] = NewVal;
                debug[3] = vario;
            }
        }
    }
}

void getAltitudePID(void)                                                     // I put this out of getEstimatedAltitude seems logical
{
    float ThrAngle;
    ThrAngle = constrain(TiltValue * 100.0f, 0, 100.0f);
    BaroP    = BaroI = BaroD = 0;                                             // Reset the Pid, create something new, or not....
    if (ThrAngle < 40 || UpsideDown) return;                                  // Don't do BaroPID if copter too tilted
    BaroP    = (int16_t)((float)cfg.P8[PIDALT] * (AltHold - EstAlt) * 0.005f);
//  BaroI    = (int16_t)((float)cfg.I8[PIDALT] * vario * 0.02f);              // That is actually a "D"
    BaroI    = (int16_t)(((float)cfg.I8[PIDALT] * vario * 0.00006f) / ACCDeltaTimeINS); // That is actually a "D"
    BaroD    = (int16_t)((float)cfg.D8[PIDALT] * (100.0f - ThrAngle) * 0.04f);// That is actually the Tiltcompensation
}
#endif

/*
/////// GPS INS TESTCODE
//  Testcode
    static uint32_t previous5HzT = 0;
		flthead = 0;                                                              // if no mag do bodyframe below
//  Testcode
    int16_t knob = constrain(rcData[AUX3]-1000,0,1000);
		float  knobbi= (float)knob * 0.001f;
		debug[0] = knobbi * 1000;
	  if (currentT > previous5HzT + 200000){
        previous5HzT = currentT;
		    VelNorth     = VelNorth * knobbi;
        VelEast      = VelEast  * knobbi;
		}
		debug[1] = VelNorth;
		debug[2] = VelEast;
/////// GPS INS TESTCODE
*/
