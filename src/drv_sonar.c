#include "board.h"
#include "mw.h"

#ifdef SONAR
/*
Currently supported Sonars: (the links are just examples, I didn't buy them there)
HC-SR04  *** Warning: HC-SR04 operates at +5V *** (google that, price around 5$ range about 2m)

DaddyWalross did an I2C "pimp" of the HC-SR04 with the idea in mind to connect several modules at once:
http://fpv-community.de/showthread.php?20988-I%B2C-Sonarplatinen

Maxbotics Sonars using PWM readout (expensive stuff, probably some people have them from APM)

ToDo?:
Other I2C Sonars: SRF02 SRF08 SRF10 SRC235 (all the same protocol)
http://www.shop.robotikhardware.de/shop/catalog/product_info.php?products_id=121

SOME INFO:

UncompDivisor:
==============
+35C 352,17 m/s = 56,79 is the divisor
+20C 343,46 m/s = 29,115 us/cm so we measure double time so: 58,23 is the divisor
-10C 325,35 m/s = 61,47 is the divisor

Error chart we measure 5823us:
+35 C  102 CM
+20 C  100 CM
-10 C   95 CM
So no Temp compensation and taking "58" as divisor seems sufficient to me.
But for the precise people, here is further info:

Temp. Compensation (taken from Maxbotics):
==========================================
(Source: http://www.maxbotix.com/documents/Temperature_Compensation.pdf)

40C to +65C (with limited operation to +85C).
Temperature Compensation that uses the time of flight in seconds, and temperature in degrees centigrade and yields the
distance in meters works for all of our products.

Dm = TOF *((20.05*SQRT(Tc+273.15))/2)
TOF is the measured Time Of Flight in seconds,
Tc is the ambient temperature in degrees C,
Dm is the distance in meters.

For 23 degrees C and 0.0058 seconds (or 5.8mS) the
distance calculates to 1.0006 meter.
If using the Serial output, first convert the distance reported by the sensor to TOF by using 147uS per inch
(TOF = inches * 1.47E-4) or 58uS per cm (TOF = cm * 5.8E-5) and then insert the TOF into the above formula.
*/

#define SONARDW_ADDRESS      0x20       // DaddyW I2C SONAR, Standard address 0x20 to 0x27  7bit!!!
#define SONARDW_DISTANCE_OUT 0x32       // NOT USED HERE #define SONARDW_ADC_OUT           0x33
#define SR04Cycle       60              // Cycletime in ms
#define MaxboticsCycle 100              // Cycletime in ms
#define DaddyWCycle     60              // Since DaddyW Sonar uses SR04 we assume that updaterate
#define UncompDivisor   58              // Temp Compensation currently not implemented

static uint16_t  trigger_pin;
static uint16_t  echo_pin;
static uint32_t  exti_line;
static uint8_t   exti_pin_source;
static IRQn_Type exti_irqn;
static uint32_t  last_measurement = 0;  // Force 1st measurement in XXX_get_distance
static uint16_t  PulseLimitInUs;
static volatile  int32_t result;
static volatile  uint32_t calltime;

void ECHO_EXTI_IRQHandler(void)
{
    static uint32_t timing_start = 0;
    uint32_t pulse_duration;
    EXTI_ClearITPendingBit(exti_line);                                  // This must be done here, otherwise gcc will optimize it to being useless, dunno why, was mentioned on a website.

    calltime = micros();
    if(digitalIn(GPIOB, echo_pin)) timing_start = calltime;             // Up flank detected
    else
    {
        pulse_duration = calltime - timing_start;
        if (pulse_duration > 174 && pulse_duration < PulseLimitInUs && timing_start) result = pulse_duration; // something like 3 cm - X (X selected by sonartype)
        else result = 0;
        timing_start = 0;
    }
}

void EXTI1_IRQHandler(void)
{
    ECHO_EXTI_IRQHandler();
}

void EXTI9_5_IRQHandler(void)
{
    ECHO_EXTI_IRQHandler();
}

bool Snr_init(void)
{
	  uint8_t          bufdaddy[2];                               // Dummy for i2c testread
    gpio_config_t    gpio;
    EXTI_InitTypeDef EXTIInit;

    // enable AFIO for EXTI support - already done is drv_system.c
    // RCC_APB2PeriphClockCmd(RCC_APB2Periph_AFIO | RCC_APB2Periph, ENABLE); 
    // cfg.snr_type = X;  0 = PWM56, 1 = RC78, 2 = I2C (DaddyWalross), 3 = MBPWM56, 4 = MBRC78
    PulseLimitInUs = 62000;                                     // Datasheet Limit 62ms for Maxbotix
    switch (cfg.snr_type)                                       // Switch according physical connection pins
    {
    case 0:                                                     // case sonar_pwm56
        PulseLimitInUs  = 24000;                                // HC-SR04 Limit 413 cm
    case 3:
        if (NumberOfMotors > 4 || feature(FEATURE_SERVO_TILT)) return false; // End here, if PWM5/6 Sonar not possible
        trigger_pin     = Pin_8;                                // PWM5 (PB8) - 5v tolerant
        echo_pin        = Pin_9;                                // PWM6 (PB9) - 5v tolerant
        exti_line       = EXTI_Line9;
        exti_pin_source = GPIO_PinSource9;
        exti_irqn       = EXTI9_5_IRQn;
        break;
    case 1:                                                     // case sonar_rc78
        PulseLimitInUs  = 24000;                                // HC-SR04 Limit 413 cm
    case 4:
        if (!feature(FEATURE_PPM)) return false;                // End here, if rc7/8 Sonar not possible, further motornumber checks missing
        trigger_pin     = Pin_0;                                // RX7 (PB0) - only 3.3v ( add a 1K Ohms resistor )
        echo_pin        = Pin_1;                                // RX8 (PB1) - only 3.3v ( add a 1K Ohms resistor )
        exti_line       = EXTI_Line1;
        exti_pin_source = GPIO_PinSource1;
        exti_irqn       = EXTI1_IRQn;
        break;
    case 2:                                                     // case sonar_i2cDW Deal with I2C daddy walross sonar
        delay(1000);                                            // sleep for 1000ms to startup sonar
        return i2cRead(SONARDW_ADDRESS, SONARDW_DISTANCE_OUT, 2, bufdaddy); // End this here
    }

    // tp - trigger pin
    gpio.pin   = trigger_pin;
    gpio.mode  = Mode_Out_PP;
    gpio.speed = Speed_2MHz;
    gpioInit(GPIOB, &gpio);
    // ep - echo pin
    gpio.pin   = echo_pin;
    gpio.mode  = Mode_IN_FLOATING;
    gpioInit(GPIOB, &gpio);
    // setup external interrupt on echo pin
    gpioExtiLineConfig(GPIO_PortSourceGPIOB, exti_pin_source);
    
    EXTI_ClearITPendingBit(exti_line);
    EXTIInit.EXTI_Line    = exti_line;
    EXTIInit.EXTI_Mode    = EXTI_Mode_Interrupt;
    EXTIInit.EXTI_Trigger = EXTI_Trigger_Rising_Falling;
    EXTIInit.EXTI_LineCmd = ENABLE;
    EXTI_Init(&EXTIInit);
    NVIC_EnableIRQ(exti_irqn);
    return true;
}

static void CheckDisconnect(void)                               // Important! Check for disconnect in flight to avoid stall value do a copter rocketjump
{
    if ((micros() - calltime) > 300000)                         // Interrupt has stalled for ca. 300ms, disconnect error.
    {
        sensorsClear(SENSOR_SONAR);
        SonarStatus = 0;
    }
}

bool hcsr04_get_distancePWM(volatile int32_t *distance)         // HC-SR04
{
    uint32_t current_time = millis();
    static bool inidone = false;
    if(current_time < (last_measurement + SR04Cycle)) return false; // repeat interval should be greater 60ms. Avoid interference between measurements.
    last_measurement = current_time;
    *distance = result / UncompDivisor;
    digitalHi(GPIOB, trigger_pin);                              // Trigger new "Ping"
    delayMicroseconds(11);                                      // The width of trig signal must be greater than 10us
    digitalLo(GPIOB, trigger_pin);
    if (inidone) CheckDisconnect();                             // Skip first trigger/read cycle for disconnect test
    else inidone = true;
    return true;                                                // Always return true when time is up, even with errorvalue of 0
}

bool MaxBotix_get_distancePWM(volatile int32_t *distance)       // Maxbotics PWM support. Tested on MB1200 XL-MaxSonar-EZ0.
{                                                               // Just Connect Echopin!!
    uint32_t current_time = millis();
    if(current_time < (last_measurement + MaxboticsCycle)) return false;// MaxBtx needs min 99 ms here
    last_measurement = current_time;
    *distance = result / UncompDivisor;
    CheckDisconnect();                                          // Check directly, because no triggering is needed and interrupt is running for seconds already (sonar disabled during groundaltinit)
    return true;                                                // Always return true when time is up, even with errorvalue of 0
}

bool DaddyW_get_i2c_distance(int32_t *distance)
{
    uint8_t  buf[2];
    int16_t  temp;
    uint32_t current_time = millis();
    if(current_time < (last_measurement + DaddyWCycle)) return false; // Hopefully DaddyW Sonar works with that timing
    last_measurement = current_time;
    if(i2cRead(SONARDW_ADDRESS, SONARDW_DISTANCE_OUT, 2, buf))
    {
        temp = (int16_t)((buf[1] << 8) | buf[0]);
        *distance = (int32_t)(temp / UncompDivisor);
	  }
    else *distance = 0;                                         // Error, Sonar not responding
    return true;                                                // Always return true when time is up, even with errorvalue of 0
}
#endif
