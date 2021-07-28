package colorlog


abstract class Presets {

    abstract String getReset()
    abstract String getDim()
    abstract String getBlack()
    abstract String getGreen()
    abstract String getYellow()
    abstract String getYellowBold()
    abstract String getBlue()
    abstract String getPurple()
    abstract String getCyan()
    abstract String getWhite()
    abstract String getRed()

    Map palette(String args) {
        return this.getProperties()
                .findAll{ key,value ->
                    key != 'class' && key != 'map' && key != 'empty' && key != 'map'}
    }
}